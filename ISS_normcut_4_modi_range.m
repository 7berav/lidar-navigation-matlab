% 분할 
% kd tree -> knn

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8,:);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 &  P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13,:);

TR = stlread('Apollo Soyuz.stl');%
PC = mesh2pc(TR, Numpoints = 90000,SamplingMethod='UniformSampling');%2022, 2024에도 없는 함수임 2025 사용
%P = TR.Points;


if exist('DATA4.mat','file')
    S_Data = load('DATA4.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA4.mat 파일이 없습니다.');
end
%P_use = [P1(1:7:end,:) ; P2(1:6:end,:);P3(1:11:end,:) ; P4(1:9:end,:)];

%{
P_use = [P3(1:11:end,:) ; P4(1:9:end,:)];

Pi = [Data_Data(5).P ; Data_Data(6).P ;  Data_Data(7).P; Data_Data(8).P; Data_Data(9).P; Data_Data(10).P];
P_use = [P_use ; Pi(1:10:end,:)];
Pi2 = [Data_Data(11).P ; Data_Data(12).P ; Data_Data(13).P ; Data_Data(21).P; Data_Data(22).P; Data_Data(23).P];
P_use = [P_use ; Pi2(1:11:end,:)];
Pi3 = [Data_Data(24).P ; Data_Data(25).P ];
P_use = [P_use ; Pi3(1:11:end,:)];
%}
%P_use = [P1(1:2:end,:) ; P2(1:2:end,:);P3(1:2:end,:) ; P4(1:2:end,:)];

P_use = PC(1:1:end,:);

%% --- (옵션) mutual-kNN 전에 반경 기반 얇게 깎기 (KD-tree, r=0.6) ---
r_thin = 0.15;                     % 반경 (필요시 조정)
N0     = size(P_use,1);

% KD-tree (거리: 유클리드)
Mdl0 = createns(P_use, 'NSMethod','kdtree', 'Distance','euclidean');

% 무작위 순서로 훑으면서 반경 r 내 이웃을 억제
ord    = randperm(N0);
keep   = false(N0,1);
blocked= false(N0,1);

for t = 1:N0
    i = ord(t);
    if blocked(i), continue; end
    keep(i) = true;

    % i와 거리 r_thin 이하인 점들을 한 번에 차단
    nbrs = rangesearch(Mdl0, P_use(i,:), r_thin);
    blocked(nbrs{1}) = true;
end

% 축소 결과 및 인덱스 매핑
idx_keep = find(keep);                  % 원본 인덱스 -> 유지되는 점 인덱스
P_use    = P_use(idx_keep, :);          % 축소된 점군
N1       = size(P_use,1);

old2new = zeros(N0,1);                  % 원본 -> 축소 인덱스
old2new(idx_keep) = 1:N1;

% (선택) 유지율 확인
fprintf('[thin] kept %d / %d (%.1f%%), r=%.3f\n', N1, N0, 100*N1/N0, r_thin);


%%

% 1) KD-tree 생성 (Statistics/ML Toolbox)
Mdl = createns(P_use, 'NSMethod','kdtree', 'Distance','euclidean');

k = 30;                                % 최근접 이웃 수
[idx, dist] = knnsearch(Mdl, P_use, 'K', k+1);
idx  = idx(:,2:end);                    % 자기 자신 제거
dist = dist(:,2:end);


% kNN 간선 리스트 (유향) → 그래프 만들 때 대칭화 예정
rows_all = repmat((1:size(P_use,1))', k, 1);  % N*k x 1
cols_all = idx(:);                % N*k x 1
di_all = dist(:);   

%mutual kNN
N = size(P_use,1);
A = sparse(rows_all, cols_all, 1, N, N);   % 유향 존재행렬
M = A & A.';                       % 상호 이웃만 유지
mask = M(sub2ind([N,N], rows_all, cols_all));
rows = rows_all(mask);            % E x 1
cols = cols_all(mask);            % E x 1
di    = di_all(mask); 

gamma = 0.90; 
sigma_i = median(dist(:));
tau = dist(:, 3);
sigma = 0.95 * ones(length(tau),1); 
tau_rows = tau(rows);
tau_cols = tau(cols);
%w_sim = exp( -(di.^2) ./ ( (tau(rows) .* tau(cols)) + eps ) ); 
%w_sim = exp( -(di.^2) ./ ( (sigma(rows) .* sigma(cols)) + eps ) ); 
w_sim = exp( -(di.^2) ./ ( (tau(rows).^gamma) .* (tau(cols).^gamma) .* (sigma_i.^(2*(1-gamma))) + eps ) ); 

% A : graph cut!  -->  Normalized Cut (spectral, K=2)
% 1) 무향 그래프 G는 위에서 생성됨
G = graph(rows, cols, w_sim, size(P_use,1));
G = simplify(G,'min');
comp = conncomp(G);
labels0 = uint32(conncomp(G)).';   % 연결 성분 ID가 곧 클러스터


W  = adjacency(G, 'weighted');    % N x N sparse
W  = max(W, W.');                  % 수치적 비대칭 보정
%%
%2-2) 개수기준 필터링
% 다만 지금 필터링 해서 아이디 바뀌는 문제가 있어서, 이거를 P_use 단에서 줄이는게 맞을지, 아니면 
% 새 점 번호 2 붙여서 만들었는데 노테이션 어떻게 하는게 좋을까
min_size   = 150;
comp_size  = accumarray(comp(:), 1);    % 크기: C×1
C = max(comp);
accumarray(comp(:), 1)
keep_comp  = comp_size >= min_size;     % true/false by component

keep_nodes = keep_comp(comp(:));   
P_use2 = P_use(keep_nodes, :);%줄임
W2 = W(keep_nodes, keep_nodes); %이것도 줄임
G2 = graph(W2, 'upper');
N2  = size(P_use2,1);
N_old    = size(P_use,1);
idx_keep = find(keep_nodes);

old2new  = zeros(N_old,1);
old2new(idx_keep) = 1:N2;  

%엣지도 마스킹 
edge_mask = keep_nodes(rows) & keep_nodes(cols);
rows2 = old2new(rows(edge_mask));      % new id 간선
cols2 = old2new(cols(edge_mask));
w2    = w_sim(edge_mask);

comp2 = conncomp(G2);
C2 = max(comp2);
accumarray(comp2(:), 1) % 개수 확인
cmap = lines(C2);

figure(11); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85);
hold on
for cid = 1:C
        idx = (comp2 == cid);
        if any(idx)
            scatter3(P_use2(idx,1), P_use2(idx,2), P_use2(idx,3), ...
            3, cmap(cid,:), 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.8);
        end
end

axis equal
%xlim([-10.5 10.5]); ylim([-15.5 15.5]); zlim([-30.5 0.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar
hold off
%%
%{
% 시각화용 서브샘플
% 여기도 새 넘버링 사용
Eidx = (1:numel(rows2))';                  % 유향 간선 인덱스
max_edges_to_draw = 30000;                % 필요시 조정
if numel(Eidx) > max_edges_to_draw
    Eidx = Eidx(randperm(numel(Eidx), max_edges_to_draw));
end

ru = rows2(Eidx); cu = cols2(Eidx);

X = [P_use2(ru,1) P_use2(cu,1) NaN(size(ru))];
Y = [P_use2(ru,2) P_use2(cu,2) NaN(size(ru))];
Z = [P_use2(ru,3) P_use2(cu,3) NaN(size(ru))];

cmap = parula(256);  % 또는 jet(256)
colormap(cmap);
C = interp1(linspace(0, 1, size(cmap,1)), cmap, w2(Eidx), 'linear', 'extrap');

figure(3);  axis equal; box on;
scatter3(P_use2(:,1), P_use2(:,2), P_use2(:,3), 1, [0.8 0.8 0.8], 'filled');hold on;
%line(X', Y', Z', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.4]);
for i = 1: numel(ru)
    line(X(i,:), Y(i,:), Z(i,:), 'Color', C(i,:), 'LineWidth', 0.3);
end
colorbar; caxis([0 1]);
view([0.3,1,0.3]); axis equal;
title(sprintf('kNN links (drawn %d / %d)', numel(ru), numel(rows)));
%}
%%
% 3) 정규화 라플라시안 L_sym = I - D^{-1/2} W D^{-1/2}
d         = full(sum(W2,2));                       % 차수벡터 d_i = Σ_j W_ij
DinvSqrt  = spdiags(1./sqrt(d + eps), 0, N2, N2);   % D^{-1/2}
Lsym      = speye(N2) - DinvSqrt * W2 * DinvSqrt;   % L_sym
Lsym = (Lsym + Lsym.')/2; 
% 4) L_sym의 가장 작은 고유값 2개의 고유벡터(U)를 계산
%    (두 번째 고유벡터 = Fiedler vector)
[U, lamb] = eigs(Lsym,45, 'smallestabs', 'Tolerance', 1e-4, 'MaxIterations', 500);


figure(10); clf;
lambda_vals = diag(lamb);              % 15개의 고유값 추출
plot(1:length(lambda_vals), lambda_vals, 'o-', 'LineWidth', 1.5);
xlabel('Eigenvalue index');
ylabel('Eigenvalue (λ)');
title('L_{sym} spectrum (smallest 15 eigenvalues)');
grid on; box on;



%%
cluster = 5; targetcut = 2;
fied = U(:,cluster+targetcut);
Y = DinvSqrt * U;

fied    = DinvSqrt * fied;                % y = D^{-1/2} z  (degree 보정)


%양수 음수로 구분
margin = 0.00053;  % 예: 0.01 ~ 0.02 사이로 조정
isPlus = fied >= 0;    % 0으로 간주되는 노드
%labels = uint32(1 + isPlus); % 1: 비0인 군, 2: 0 군
labels = zeros(size(fied), 'uint32');
labels(fied < -margin)        = 1;  % 1: 음수 쪽
labels(abs(fied) <= margin)   = 2;  % 2: 0 근처(마진 밴드)
labels(fied >  margin)        = 3;




% 6) 셀 배열로 추출 (clusters{1}, clusters{2})
K = 3;
clusters = cell(1,K);
for kcid = 1:K
    clusters{kcid} = P_use2(labels==kcid, :);
end
clusters = clusters(~cellfun('isempty',clusters));
cmap = lines(K);
[usorted, order] = sort(fied);
labels_sorted = labels(order);

figure(9); clf; hold on;
for kcid = 1:K
    idx = (labels_sorted == kcid);
    scatter(find(idx), usorted(idx), 15, cmap(kcid,:), 'filled');
end

xlabel('Node rank (sorted by Y(:,2))');
ylabel('U(:,2) value');
title('Sorted Fiedler vector');
grid on; box on;
hold off;


figure(12); clf;
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; 
if K > 0
    cmap = lines(K);
    for kcid = 1:K
        idx = (labels == kcid);
        if any(idx)
            scatter3(P_use2(idx,1), P_use2(idx,2), P_use2(idx,3), ...
                     3, cmap(kcid,:), 'filled', ...
                     'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4);
        end
    end
end

hold off
axis equal
%xlim([-10.5 10.5]); ylim([-15.5 15.5]); zlim([-30.5 0.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');

pq   = [cluster+2 cluster+4];                 % <- 여기만 바꾸세요. [3 5]로 하면 3:5
idx3 = pq(1):pq(2);

X3 = Y(:, idx3);              % 선택한 3개 고유벡터
% (권장) 행 정규화: 각 노드를 단위벡터로
%rn = sqrt(sum(X3.^2,2)) + eps;
%X3 = X3 ./ rn;

figure(15); clf; hold on;
if exist('labels','var') && numel(labels)==size(U,1) && numel(unique(labels))>1
    scatter3(X3(:,1), X3(:,2), X3(:,3), 4, double(labels), 'filled');
    colormap(lines(numel(unique(labels))));
else
    scatter3(X3(:,1), X3(:,2), X3(:,3), 2, 'filled');
end
xlabel(sprintf('u_{%d}', idx3(1))); ylabel(sprintf('u_{%d}', idx3(2))); zlabel(sprintf('u_{%d}', idx3(3)));

title(sprintf('Spectral embedding: eigenvectors %d:%d', idx3(1), idx3(3)));
grid on; axis equal; view([1,-1,1]); hold off;
xlim([-0.02 0.02]); ylim([-0.02 0.02]); zlim([-0.02 0.02]);   % 필요에 맞게

%% make new cluster
ndim  = cluster + 15;               % 사용할 고유벡터(축) 개수
Y_use     = Y(:, cluster+1:ndim);                      % 앞에서부터 ndim축 사용 (필요시 조정)
%벡터스페이스 mst 만들기
Ny  = size(Y_use, 1);
ky  = 40;    % 필요시 조정
% 1) kNN (자기 자신 제외)

Mdly = createns(Y_use, 'NSMethod','kdtree', 'Distance','euclidean');

[idxNN_y, distNN_y] = knnsearch(Mdly, Y_use, 'K', ky+1);
idxNN_y  = idxNN_y(:, 2:end);          % Ny x ky
distNN_y = distNN_y(:, 2:end);         % Ny x ky


% 유향 kNN → 희소행렬 
I_y  = repmat((1:Ny).', ky, 1);
J_y  = idxNN_y(:);
W_y  = distNN_y(:);                    % 거리(= MST 가중치)
S_y  = sparse(I_y, J_y, W_y, Ny, Ny);  % 유향: 거리 가중치
S_yu = max(S_y, S_y.');
% 3) 무향 그래프 및 MST
G_knn_y = graph(S_yu, 'upper');

compGknn = conncomp(G_knn_y);                     % 1..m
m = max(compGknn);
sg = accumarray(conncomp(G_knn_y).',1);

T_cells = cell(m,1);
for i = 1:m
    Vi = find(compGknn == i);
    Gi = subgraph(G_knn_y, Vi);
    Ti = minspantree(Gi, 'Method','sparse');   % 성분 i의 MST
    % 원래 인덱스로 되돌리기 위해 노드 인덱스 보관
    Ti.Nodes.orig = Vi(:);                     % 메모: plot에는 불필요하지만 편의상 보관
    T_cells{i} = Ti;
end

%% cut
K_target_y  = 80;                                   % 목표 군집 수
N = Ny;  % 미정의였던 N을 Ny로 고정
need_cut = max(0, K_target_y - m);
q_thr = 0.96; % 상위 1% 기준치

%{
Wv = W_all(isfinite(W_all));

% 분위수
q95 = quantile(Wv, 0.95);
q99 = quantile(Wv, 0.99);

% 히스토그램 (확률밀도 PDF 기준)
figure(910); clf; hold on;
histogram(Wv, 'NumBins', 80, 'Normalization', 'pdf', ...
          'FaceColor', [0.65 0.75 0.95], 'EdgeColor', 'none');

% 분위선 그리기
xline(q95, 'r--', 'LineWidth', 1.3, 'Label', '95%', 'LabelVerticalAlignment','bottom');
xline(q99, 'm-.', 'LineWidth', 1.3, 'Label', '99%', 'LabelVerticalAlignment','bottom');

% 간단 표기(선택)
yl = ylim;
text(q95, yl(2)*0.92, sprintf('%.4g', q95), 'Color','r', 'HorizontalAlignment','left');
text(q99, yl(2)*0.80, sprintf('%.4g', q99), 'Color','m', 'HorizontalAlignment','left');

xlabel('Edge weight (W)');
ylabel('PDF');
title('Histogram of edge weights with 95% / 99% quantiles');
grid on; box on; hold off;
%}

% 1) 모든 성분의 간선(전역 인덱스, 거리=Weight) 모으기
all_edges = [];
for i = 1:m
    Ti = T_cells{i};
    E  = Ti.Edges.EndNodes;     % 성분 내부 인덱스
    W  = Ti.Edges.Weight;       % 간선 길이(거리)
    Vi = Ti.Nodes.orig;         % 원래 노드 인덱스
    all_edges = [all_edges; [Vi(E(:,1)), Vi(E(:,2)), W]];
end
W_all = all_edges(:,3);

% 2) 분위수 임계 계산 & 컷 대상 간선 선정
tau = quantile(W_all, q_thr)        % 임계 길이
cut_list = all_edges(W_all >= tau, 1:2);

% 3) 성분별로 해당 간선 제거
if ~isempty(cut_list)
    for i = 1:m
        Ti = T_cells{i};
        Vi = Ti.Nodes.orig;
        mask_i = ismember(cut_list(:,1), Vi) & ismember(cut_list(:,2), Vi);
        if any(mask_i)
            map = zeros(N,1); map(Vi) = 1:numel(Vi);     % 전역→성분 내부 인덱스
            Ecut_i = [map(cut_list(mask_i,1)), map(cut_list(mask_i,2))];
            Ti = rmedge(Ti, Ecut_i(:,1), Ecut_i(:,2));
        end
        T_cells{i} = Ti;
    end
end
labels_all = zeros(N,1,'uint32');
offset = 0;
for i = 1:m
    Ti = T_cells{i};
    Vi = Ti.Nodes.orig;
    li = conncomp(Ti).';
    labels_all(Vi) = uint32(li + offset);
    offset = offset + max(li);
end
%최소 사이즈 기준
min_size = 80;              
K_now = double(max(labels_all));
sz = accumarray(double(labels_all(labels_all>0)), 1, [K_now, 1]);
keep = find(sz >= min_size);
labels_display = zeros(N,1,'uint32');
if ~isempty(keep)
    % 연속 라벨로 재맵핑
    remap = zeros(K_now,1,'uint32'); remap(keep) = 1:numel(keep);
    labels_display(labels_all>0) = remap(labels_all(labels_all>0));
end
sz_cut = accumarray(double(labels_display(labels_display>0)), 1);
%%


pq   = [3 5];                 % <- 여기만 바꾸세요. [3 5]로 하면 3:5
idx3 = pq(1):pq(2);

N = Ny;  % 미정의였던 N을 Ny로 고정
E_global = [];
W_global = [];
for i = 1:m
    Ti = T_cells{i};
    Vi = Ti.Nodes.orig;              % 원래 전역 인덱스
    E  = Ti.Edges.EndNodes;          % 성분 내부 인덱스
    W  = Ti.Edges.Weight;            % 가중치(거리)
    E_global = [E_global; [Vi(E(:,1)), Vi(E(:,2))]];
    W_global = [W_global; W(:)];
end

% 포레스트 그래프 (전역 노드 수 N)
% (가중치 보존이 필요하면 세 번째 인자로 W_global을 넣습니다)
T_all = graph(E_global(:,1), E_global(:,2), W_global, N);

d = size(Y_use, 2);
imax = max(idx3);
if imax > d
    Y_plot = [Y_use, zeros(size(Y_use,1), imax - d)];
else
    Y_plot = Y_use(:, idx3);
end

figure(200); clf; hold on;
% 그래프 plot (노드 좌표는 Y_plot, 간선은 포레스트 간선)
h = plot(T_all, ...
         'XData', Y_plot(:,1), ...
         'YData', Y_plot(:,2), ...
         'ZData', Y_plot(:,3), ...
         'LineWidth', 0.6, ...
         'MarkerSize', 0.1, ... 
         'NodeLabel', {});
view(3); axis equal; grid on;
% 노드 색: degree
deg = degree(T_all);
h.NodeCData = deg;
deg_min = min(deg);
deg_max = max(deg);
eps_n = 1e-12;                           % 분모 0 방지용

sz_min = 2;                              
sz_max = 8;
w = (deg - deg_min) ./ max(deg_max - deg_min, eps_n);
sz = sz_min + (sz_max - sz_min) .* w;

sc = scatter3(Y_plot(:,1), Y_plot(:,2), Y_plot(:,3), ...
              sz, deg, 'filled', ...
              'MarkerEdgeColor', 'none','MarkerFaceAlpha',0.4, 'LineWidth', 0.2);

colormap(jet); colorbar;
colormap(jet); colorbar;
%title('MSF on Y\_use (node color = degree)');
xlabel(sprintf('y_{%d}', idx3(1))); ylabel(sprintf('y_{%d}', idx3(2))); zlabel(sprintf('y_{%d}', idx3(3)));
xlim([-0.02 0.02]); ylim([-0.02 0.02]); zlim([-0.02 0.02]); 
view([1 1 1]);
hold off;

%% eigenvector space plot

% 전제: Y_plot은 이미 idx3 축을 반영해 구성되어 있음 (figure 200에서 만든 동일 변수 사용)
figure(16); clf; hold on;
K_disp = double(max(labels_display));
cmap = lines(K_disp);
if exist('labels_display','var') && any(labels_display>0) && numel(unique(labels_display(labels_display>0)))>0
    scatter3(Y_plot(:,1), Y_plot(:,2), Y_plot(:,3), ...
             3, double(labels_display), 'filled', ...
             'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.3);
    colormap(cmap);
else
    scatter3(Y_plot(:,1), Y_plot(:,2), Y_plot(:,3), ...
             2, [0.3 0.3 0.3], 'filled', ...
             'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
end
xlabel(sprintf('u_{%d}', idx3(1)));
ylabel(sprintf('u_{%d}', idx3(2))); 
zlabel(sprintf('u_{%d}', idx3(3))); 
%title(sprintf('Spectral embedding (dims [%s])', num2str(idx3)));
grid on; axis equal; view([1 1 1]);
xlim([-0.02 0.02]); ylim([-0.02 0.02]); zlim([-0.02 0.02]); 

hold off;

%% mst draw

idx_keep = [];
if exist('labels_display','var'), idx_keep = find(labels_display>0); end

figure(300); clf; axis equal; box on; hold on;
if ~isempty(idx_keep)
    % 유지 노드 서브그래프
    T_sub = subgraph(T_all, idx_keep);

    % 좌표: 반드시 Y_plot(=idx3 축 반영) 사용
    Xc = Y_plot(idx_keep,1); 
    Yc = Y_plot(idx_keep,2); 
    Zc = Y_plot(idx_keep,3);

    % 배경 노드(유지분) 산점도: 작게 + 반투명
    scatter3(Xc, Yc, Zc, 1, [0.6 0.6 0.6], 'filled', ...
             'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.25);

    if numedges(T_sub) > 0
        eu = T_sub.Edges.EndNodes(:,1);
        ev = T_sub.Edges.EndNodes(:,2);
        Wk = T_sub.Edges.Weight;
        w01 = (Wk - min(Wk)) ./ max(eps, max(Wk)-min(Wk));   % [0,1]

        cm = parula(256);
        for t = 1:numel(Wk)
            c = cm( 1 + floor( w01(t)*(size(cm,1)-1) ), : );
            line([Xc(eu(t)) Xc(ev(t))], ...
                 [Yc(eu(t)) Yc(ev(t))], ...
                 [Zc(eu(t)) Zc(ev(t))], ...
                 'Color', c, 'LineWidth', 0.9);
        end
        colorbar;  % 간선 색 스케일 확인용
    end
    title(sprintf('MST (kept nodes only), dims [%s]', num2str(idx3))); view([1 1 1]);
else
    text(0.1,0.5,'labels\_display > 0 노드가 없습니다.','FontSize',12);
    axis off;
end
view([1 1 1]);
%xlim([-0.010 -0.000]);ylim([-0.005 +0.006]);zlim([-0.005 0.005]);
xlim([-0.015 0.015]); ylim([-0.015 0.015]); zlim([-0.015 0.015]); 

grid on; hold off;

%% 분할 그리기 (figure 12) — labels_display 기준
figure(13); clf; hold on;

% 0라벨(버린 군집) 배경
idx_drop = (labels_display == 0);
if any(idx_drop)
    scatter3(P_use(idx_drop,1), P_use(idx_drop,2), P_use(idx_drop,3), ...
             1, [0.55 0.45 0.45], 'filled', ...
             'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);
end

% 유지된 군집(1..K_disp) 색칠
for k = 1:K_disp
    idx = (labels_display == k);
    if any(idx)
        scatter3(P_use(idx,1), P_use(idx,2), P_use(idx,3), ...
                 3, cmap(k,:), 'filled', ...
                 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.6);
    end
end

hold off; grid on; view([1 1 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-30.5 30.5]); ylim([-25.5 25.5]); zlim([-30.5 20.5]);   % 필요에 맞게
axis equal; 
title('Clusters on original point cloud (after min\_size)');
%% saving into mat 저장코드

Kmix = K_disp;

% 3D 점 좌표 선택: P_use(:,1:3) 을 기본으로 사용
P_xyz = P_use(:,1:3);

% --- Data 구조체 및 점군 셀 배열 생성 ---
Data   = repmat(struct('P',[],'Disp',[],'Beta',[],'InlierMask',[]), 1, K_disp);
P_cell  = cell(1, K_disp);
counts  = zeros(1, K_disp);

for k = 1:K_disp
    idx_k = (double(labels_display) == k);
    Pk    = P_xyz(idx_k, :);             % N_k × 3
    counts(k) = size(Pk,1);

    Data(k).P          = Pk;
    Data(k).Disp       = [];            % RANSAC에서 채우십시오 (1×3 또는 3×1)
    Data(k).Beta       = [];            % RANSAC에서 채우십시오
    Data(k).InlierMask = [];            % RANSAC에서 채우십시오

    P_cell{k} = Pk;                      % 점군만 저장하는 셀
end

% 간단 보고
fprintf('[Data8 저장] 군집 수 K = %d, 군집별 점 개수(min/median/max) = %d / %d / %d\n', ...
        K_disp, min(counts), median(counts), max(counts));

% --- 저장 ---
save('image_cut/DATA8.mat', 'Data', 'P_cell', 'labels', 'Kmix', '-v7.3');
fprintf('DATA8.mat 저장 완료: 변수 {Data6, P_cell, labels, Kmix}\n');
