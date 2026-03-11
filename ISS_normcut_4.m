% 분할 
% kd tree -> knn

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8,:);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 &  P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13,:);


if exist('DATA4.mat','file')
    S_Data = load('DATA4.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA4.mat 파일이 없습니다.');
end
P_use = [P1(1:10:end,:) ; P2(1:9:end,:);P3(1:7:end,:) ; P4(1:10:end,:)];


Pi = [Data_Data(5).P ; Data_Data(6).P ;  Data_Data(7).P];
P_use = [P_use ; Pi(1:7:end,:)];
%P_use = [P1(1:2:end,:) ; P2(1:2:end,:);P3(1:2:end,:) ; P4(1:2:end,:)];
%P_use = P(1:10:end,:);
% 1) KD-tree 생성 (Statistics/ML Toolbox)
Mdl = createns(P_use, 'NSMethod','kdtree', 'Distance','euclidean');

k = 19;                                % 최근접 이웃 수
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


sigma_i = median(dist(:));
tau = dist(:, 3);
tau_rows = tau(rows);
tau_cols = tau(cols);
%w_sim = exp(-(di.^2) ./ ((tau(rows_all) .* tau(cols_all))  + eps));   % 유사도 in (0,1]
w_sim = exp( -(di.^2) ./ ( (tau(rows) .* tau(cols)) + eps ) ); 

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
C = max(comp2)
accumarray(comp2(:), 1)
cmap = lines(C);

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

%%
% 3) 정규화 라플라시안 L_sym = I - D^{-1/2} W D^{-1/2}
d         = full(sum(W2,2));                       % 차수벡터 d_i = Σ_j W_ij
DinvSqrt  = spdiags(1./sqrt(d + eps), 0, N2, N2);   % D^{-1/2}
Lsym      = speye(N2) - DinvSqrt * W2 * DinvSqrt;   % L_sym
Lsym = (Lsym + Lsym.')/2; 
% 4) L_sym의 가장 작은 고유값 2개의 고유벡터(U)를 계산
%    (두 번째 고유벡터 = Fiedler vector)
[U, lamb] = eigs(Lsym,45, 'smallestabs', 'Tolerance', 1e-4, 'MaxIterations', 500);


figure(14); clf;
lambda_vals = diag(lamb);              % 15개의 고유값 추출
plot(1:length(lambda_vals), lambda_vals, 'o-', 'LineWidth', 1.5);
xlabel('Eigenvalue index');
ylabel('Eigenvalue (λ)');
title('L_{sym} spectrum (smallest 15 eigenvalues)');
grid on; box on;

%{
figure(24); clf;

lambda_vals = diag(lamb);  % lamb: 고유값이 대각선에 있는 행렬
plot(2:40, lambda_vals(2:40).' ./ ((2:40)-1), 'o-')
xlabel('k');
ylabel('\lambda_k / k');
title('L_{sym}: \lambda_k/k (k=2..7)');
ylim([0 3.0e-3])
grid on; box on; 
%}

%%
cluster = 6; targetcut = 2 ;
fied = U(:,cluster+targetcut);
Y = DinvSqrt * U;

fied    = DinvSqrt * fied;                % y = D^{-1/2} z  (degree 보정)


%양수 음수로 구분
margin = 0.0005;  % 예: 0.01 ~ 0.02 사이로 조정
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
%고유벡터 값 나열 
%{
[usorted, order] = sort(fied);
labels_sorted = labels(order);

figure(13); clf; hold on;
for kcid = 1:K
    idx = (labels_sorted == kcid);
    scatter(find(idx), usorted(idx), 15, cmap(kcid,:), 'filled');
end
xlabel('Node rank (sorted by Y(:,2))');
ylabel('U(:,2) value');
title('Sorted Fiedler vector (rank vs value, color = label)');
grid on; box on;
hold off;
%}
%% pointcloud 그리기

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
                     'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.8);
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
ndim  = cluster + 6;               % 사용할 고유벡터(축) 개수
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
sg = accumarray(conncomp(G_knn_y).',1);
% 무향 그래프(상삼각 사용)
T_mst_y = minspantree(G_knn_y,  'Method', 'sparse');
c0 = conncomp(T_mst_y).';    
C0 = max(c0);
st = accumarray(c0,1);
K_target_y  = 10;                                   % 목표 군집 수
E_uv_y  = T_mst_y.Edges.EndNodes;               % [u v] (Nx2 double)
W_e_y   = T_mst_y.Edges.Weight;                 % 간선 길이
[~, ord] = sort(W_e_y, 'descend');         % 긴 간선부터
cut_m = max(0, K_target_y-1);
if cut_m > 0
    cut_ids = ord(1:min(cut_m, numel(ord)));
    T_cut = rmedge(T_mst_y, E_uv_y(cut_ids,1), E_uv_y(cut_ids,2));
else
    T_cut = T_mst_y;
end
labels_mst = conncomp(T_cut).';  % 1..Kmix

min_size = 2;                 % 유지할 최소 클러스터 크기(변수)
K = K_target_y;                      % 혼합 성분 수
% 1) 각 라벨의 크기 계산 (labels ∈ {1..K})
sz = accumarray(double(labels_mst(:)), 1);  % K×1
keep_label_mask = sz >= min_size;         % K×1 logical
keep_labels     = find(keep_label_mask);
% 최소 사이즈 라벨링 ( 겸 다른 그룹 conncomp 결과 제외)
keep_nodes_mask   = ismember(labels_mst, keep_labels);   % N2×1 logical
labels_mst_kept   = labels_mst;                          % N2×1
labels_mst_kept(~keep_nodes_mask) = 0;                   % 버린 라벨 → 0

sz_kept = accumarray(double(categorical(labels_mst(keep_nodes_mask), keep_labels)), 1);

% (선택) 유지 노드만으로 서브그래프/라벨 재계산 (권장)
idx_keep = find(keep_nodes_mask);
T_keep   = subgraph(T_cut, idx_keep);
labels_keep_sub = conncomp(T_keep).';    
labels_mst_kept_sub = zeros(Ny,1,'uint32');
labels_mst_kept_sub(idx_keep) = uint32(labels_keep_sub);
%% test knn

%distNN_y2 = norm(Y2 )%그냥 기존 idx 쌍 사용해서 거리 구하기 
%그거로 가중치만 바꾼 G_knn_y2 만들기
% 그거로 minspantree 돌리기

%앞이랑 뒤랑 비교해서 그림그리기? cut 이전에 원래 데이터로 그림그리기/ 
Y_use;               % Ny x d (d>=1)
d   = size(Y_use,2);
if d < 2, Y_use(:,2) = 0; end  % 1D -> 2D로 보정

figure(200); clf; hold on;
if d >= 3
    % 3D plot using graph plot (supports ZData)
    h = plot(T_mst_y, 'XData', Y_use(:,2), 'YData', Y_use(:,3), 'ZData',  Y_use(:,4), ...
             'LineWidth', 0.2, 'NodeLabel', {});
    view(3);

end

% 노드 색깔을 degree로 칠하기 (선택)
deg = degree(T_mst_y);
h.NodeCData = deg;
colormap(jet); colorbar;
title('MST on Y\_use\_y (node color = degree)');
axis equal; grid on;



%% eigenvector space plot

pq   = [2 4];                 % <- 여기만 바꾸세요. [3 5]로 하면 3:5
idx3 = pq(1):pq(2);

X3 = Y_use(:, idx3); 
figure(16); clf; hold on;
if exist('labels','var') && numel(labels_mst_kept)==size(U,1) && numel(unique(labels_mst_kept))>1
    scatter3(X3(:,1), X3(:,2), X3(:,3), 4, double(labels_mst_kept), 'filled');
    colormap(lines(double(numel(unique(labels_mst_kept)))));
else
    scatter3(X3(:,1), X3(:,2), X3(:,3), 2, 'filled');
end
xlabel(sprintf('u_{%d}', idx3(1)));
ylabel(sprintf('u_{%d}', idx3(2)));
zlabel(sprintf('u_{%d}', idx3(3)));
title(sprintf('Spectral embedding (Y): eigenvectors %d:%d', idx3(1), idx3(3)));
grid on; axis equal ; view([1 1 1]); hold off;
%% mst draw
idx3 = [2 3 4];
eu = T_keep.Edges.EndNodes(:,1);
ev = T_keep.Edges.EndNodes(:,2);
Wk = T_keep.Edges.Weight;
w01 = (Wk - min(Wk)) ./ max(eps, max(Wk)-min(Wk));   % [0,1]

Xc = Y_use(idx_keep, idx3(1));
Yc = Y_use(idx_keep, idx3(2));
Zc = Y_use(idx_keep, idx3(3));

figure(300); clf; axis equal; box on; hold on;
scatter3(Xc(:), Yc(:), Zc(:),1, [0.7 0.7 0.7], 'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

cm = parula(256);
for t = 1:numel(Wk)
    c = cm( 1 + floor( w01(t)*(size(cm,1)-1) ), : );
    line([Xc(eu(t)) Xc(ev(t))], [Yc(eu(t)) Yc(ev(t))], [Zc(eu(t)) Zc(ev(t))], ...
         'Color', c, 'LineWidth', 0.5);
end
title('MST (kept nodes only)'); view([1 1 1]); colorbar; hold off;
%% 분할 그리기

figure(12); clf;
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; 
if K_target_y > 0
    cmap = lines(K_target_y);
    for kcid = 1:K_target_y
        if ~keep_label_mask(kcid), continue; end 
        idx = (labels == kcid);
        if any(idx)
            scatter3(P_use2(idx,1), P_use2(idx,2), P_use2(idx,3), ...
                     3, cmap(kcid,:), 'filled', ...
                     'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.8);
        end
    end
end
hold off
axis equal
%xlim([-10.5 10.5]); ylim([-15.5 15.5]); zlim([-30.5 0.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');

%% saving into mat 저장코드
%{
% Kmix 결정 (gm이 있으면 우선 사용, 없으면 라벨에서 추정)
if exist('gm','var') && ~isempty(gm) && isprop(gm,'NumComponents')
    Kmix = gm.NumComponents;
else
    Kmix = double(max(labels));
end

% 3D 점 좌표 선택: P_use(:,1:3) 을 기본으로 사용
P_xyz = P_use(:,1:3);

% --- Data6 구조체 및 점군 셀 배열 생성 ---
Data   = repmat(struct('P',[],'Disp',[],'Beta',[],'InlierMask',[]), 1, Kmix);
P_cell  = cell(1, Kmix);
counts  = zeros(1, Kmix);

for k = 1:Kmix
    idx_k = (double(labels) == k);
    Pk    = P_xyz(idx_k, :);             % N_k × 3
    counts(k) = size(Pk,1);

    Data(k).P          = Pk;
    Data(k).Disp       = [];            % RANSAC에서 채우십시오 (1×3 또는 3×1)
    Data(k).Beta       = [];            % RANSAC에서 채우십시오
    Data(k).InlierMask = [];            % RANSAC에서 채우십시오

    P_cell{k} = Pk;                      % 점군만 저장하는 셀
end

% 간단 보고
fprintf('[Data7 저장] 군집 수 K = %d, 군집별 점 개수(min/median/max) = %d / %d / %d\n', ...
        Kmix, min(counts), median(counts), max(counts));

% --- 저장 ---
save('DATA7.mat', 'Data', 'P_cell', 'labels', 'Kmix', '-v7.3');
fprintf('DATA7.mat 저장 완료: 변수 {Data6, P_cell, labels, Kmix}\n');
%}