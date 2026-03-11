% 분할 
% kd tree -> knn

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8,:);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 &  P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13,:);


%P_use = [P1(1:100:end,:) ; P2(1:100:end,:);P3(1:100:end,:) ; P4(1:100:end,:)];
%P_use = [Pi(1:20:end,:)];
P_use = [P1(1:10:end,:) ; P2(1:10:end,:);P3(1:10:end,:) ; P4(1:10:end,:);P1(1:20:end,:)];
%P_use = P(1:2:end,:);
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
%sig_rows = sigma_i;
%sig_cols = sigma_i;
w_sim = exp(-(di.^2) ./ (2*(sigma_i.*sigma_i) + eps));   % 유사도 in (0,1]


%% A : graph cut!  -->  Normalized Cut (spectral, K=2)
% 1) 무향 그래프 G는 위에서 생성됨
G = graph(rows, cols, w_sim, size(P_use,1));
G = simplify(G,'min');
labels0 = uint32(conncomp(G)).';   % 연결 성분 ID가 곧 클러스터
disp(numel(unique(labels0)));
% 2) 가중 인접행렬 W (sparse) 및 대칭화
N  = size(P_use,1);
W  = adjacency(G, 'weighted');    % N x N sparse
W  = max(W, W.');                  % 수치적 비대칭 보정

%2-2) 개수기준 필터링
G0   = graph(M | M.', 'upper');        % 무향 이진 그래프
comp = conncomp(G0);                    % 1..C (각 노드의 성분 ID)
disp(numel(unique(comp)));
min_size   = 50;
comp_size  = accumarray(comp(:), 1);    % 크기: C×1
keep_comp  = comp_size >= min_size;     % true/false by component
keep_nodes = keep_comp(comp(:));   
P_use = P_use(keep_nodes, :);
W = W(keep_nodes, keep_nodes);
G = graph(W, 'upper');
N  = size(P_use,1);
%comp2 = conncomp(G);
%accumarray(comp2(:), 1);
% 3) 정규화 라플라시안 L_sym = I - D^{-1/2} W D^{-1/2}
d         = full(sum(W,2));                       % 차수벡터 d_i = Σ_j W_ij
DinvSqrt  = spdiags(1./sqrt(d + eps), 0, N, N);   % D^{-1/2}
Lsym      = speye(N) - DinvSqrt * W * DinvSqrt;   % L_sym
Lsym = (Lsym + Lsym.')/2; 
% 4) L_sym의 가장 작은 고유값 2개의 고유벡터(U)를 계산
%    (두 번째 고유벡터 = Fiedler vector)
[U, lamb] = eigs(Lsym, 14, 'smallestabs', 'Tolerance', 1e-4, 'MaxIterations', 500);


figure(14); clf;
lambda_vals = diag(lamb);              % 15개의 고유값 추출
plot(1:length(lambda_vals), lambda_vals, 'o-', 'LineWidth', 1.5);
xlabel('Eigenvalue index');
ylabel('Eigenvalue (λ)');
title('L_{sym} spectrum (smallest 15 eigenvalues)');
grid on; box on;

     
%%
cluster = 4; targetcut = 5;
fied = U(:,cluster+targetcut);
Y = DinvSqrt * U;

fied    = DinvSqrt * fied;                % y = D^{-1/2} z  (degree 보정)


% 5) 2-분할: Fiedler 벡터에 대해 k-means (K=2)
% 0과 아닌 애들로 분류 (0 고유벡터)
%tol    = 1e-8;               % 허용오차(데이터 스케일에 맞춰 조정)
%isZero = abs(fied) <= tol;    % 0으로 간주되는 노드
%labels = uint32(1 + isZero); % 1: 비0인 군, 2: 0 군

% 유사한 두 그룹으로 분류
%rng(0);
%labels = uint32(kmeans(fied, 2, 'Replicates', 5, 'MaxIter', 200));

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
    clusters{kcid} = P_use(labels==kcid, :);
end
clusters = clusters(~cellfun('isempty',clusters));
cmap = lines(2*K);
%고유벡터 값 나열 

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

%% pointcloud 그리기

figure(11); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, fied(:), 'filled', ...
         'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85);
axis equal
%xlim([-10.5 10.5]); ylim([-15.5 15.5]); zlim([-30.5 0.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar

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
            scatter3(P_use(idx,1), P_use(idx,2), P_use(idx,3), ...
                     3, cmap(kcid,:), 'filled', ...
                     'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.8);
        end
    end
end
hold off
axis equal
%xlim([-10.5 10.5]); ylim([-15.5 15.5]); zlim([-30.5 0.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');

pq   = [8 10];                 % <- 여기만 바꾸세요. [3 5]로 하면 3:5
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
grid on; axis vis3d; view(45,25); hold off;

%%
ndim  = cluster + 5;               % 사용할 고유벡터(축) 개수
Y_use     = Y(:, 1:ndim);                      % 앞에서부터 ndim축 사용 (필요시 조정)

% ---- GMM 클러스터링 ----
Kmix  = 10;                                 % 혼합 성분 수 (분리하고 싶은 군집 수)
opts  = statset('MaxIter', 330,'TolFun', 1e-5);
tic
gm = fitgmdist(Y_use, Kmix, ...
    'SharedCovariance', true, ...          % 팔 구조에 유리(편향 완화)
    'RegularizationValue', 1e-6, ...
    'Start', 'plus', ...
    'Options', opts);
toc
% 후처리: 라벨 추정
P = posterior(gm, Y_use);                      % N x Kmix
[~, lab] = max(P, [], 2);
labels = uint32(lab);  

min_size = 80;                 % 유지할 최소 클러스터 크기(변수)
K = Kmix;                      % 혼합 성분 수
% 1) 각 라벨의 크기 계산 (labels ∈ {1..K})
sz = accumarray(double(labels), 1, [K 1], @sum, 0);  % K×1
keep_label_mask = sz >= min_size;         % K×1 logical
% 3) 노드(포인트) 마스크: 자신의 라벨이 keep_label_mask에서 true인 것만 유지
%keep_nodes = keep_label_mask(double(labels));   % N×1 logical


%% eigenvector space plot

pq   = [7 9];                 % <- 여기만 바꾸세요. [3 5]로 하면 3:5
idx3 = pq(1):pq(2);

X3 = Y(:, idx3); 
figure(16); clf; hold on;
if exist('labels','var') && numel(labels)==size(U,1) && numel(unique(labels))>1
    scatter3(X3(:,1), X3(:,2), X3(:,3), 4, double(labels), 'filled');
    colormap(lines(double(numel(unique(labels)))));
else
    scatter3(X3(:,1), X3(:,2), X3(:,3), 2, 'filled');
end
xlabel(sprintf('u_{%d}', idx3(1)));
ylabel(sprintf('u_{%d}', idx3(2)));
zlabel(sprintf('u_{%d}', idx3(3)));
title(sprintf('Spectral embedding (Y): eigenvectors %d:%d', idx3(1), idx3(3)));
grid on; axis vis3d; view(45,25); hold off;

%% 분할 그리기

figure(12); clf;
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; 
if Kmix > 0
    cmap = lines(2*Kmix);
    for kcid = 1:Kmix
        if ~keep_label_mask(kcid), continue; end 
        idx = (labels == kcid);
        if any(idx)
            scatter3(P_use(idx,1), P_use(idx,2), P_use(idx,3), ...
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