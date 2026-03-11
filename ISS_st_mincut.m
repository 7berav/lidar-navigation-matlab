% 분할 
% kd tree -> knn
% gomory hu , mst 
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
i = 13;                                 % 예: 5번째 점군을 보고 싶을 때
Pi = Data_Data(i).P;  


%P_use = [P1 ; P2; P3; P4];
P_use = [Pi(1:200:50000,:)];
% 1) KD-tree 생성 (Statistics/ML Toolbox)
Mdl = createns(P_use, 'NSMethod','kdtree', 'Distance','euclidean');

k = 12;                                % 최근접 이웃 수
[idx, dist] = knnsearch(Mdl, P_use, 'K', k+1);
idx  = idx(:,2:end);                    % 자기 자신 제거
dist = dist(:,2:end);

% kNN 간선 리스트 (유향) → 그래프 만들 때 대칭화 예정
%rows = repmat((1:size(P_use,1))', k,1);
%cols = idx(:);
%d    = dist(:); 
rows_all = repmat((1:size(P_use,1))', k, 1);  % N*k x 1
cols_all = idx(:);                % N*k x 1
d_all = dist(:);   
N = size(P_use,1);
A = sparse(rows_all, cols_all, 1, N, N);   % 유향 존재행렬
M = A & A.';                       % 상호 이웃만 유지
mask = M(sub2ind([N,N], rows_all, cols_all));
rows = rows_all(mask);            % E x 1
cols = cols_all(mask);            % E x 1
d    = di_all(mask); 

sigma_i = 0.6;%median(dist, 2);
tau = dist(:, 3);
tau_rows = tau(rows);
tau_cols = tau(cols);
sig_rows = sigma_i;
sig_cols = sigma_i;
w_sim = exp(-(d.^2) ./ ((tau(rows) .* tau(cols)) + eps));   % 유사도 in (0,1]
%w_sim = exp(-(d.^2) ./ ((2*sig_rows .* sig_cols) + eps));   % 유사도 in (0,1]

%% A : graph cut!  -->  s–t min-cut (가장 먼 점 쌍을 단말로)
% 1) 무향 그래프(가중: 거리)
G = graph(rows, cols, w_sim, size(P_use,1));
G = simplify(G, 'min');                 % 중복 간선 최소 가중치만 유지

% 2) 그래프 최단거리 기준으로 "대략" 가장 먼 두 노드 s,t 선택
%    (정확도 중요치 않다고 하셨으니 전체 distances로 간단히)
%[~, s] = min(P_use(:,3));   % z 최솟값 노드
%[~, t] = max(P_use(:,3));   % z 최댓값 노드
q1 = prctile(P_use(:,3), 25);
q3 = prctile(P_use(:,3), 75);

[~,s] = min(abs(P_use(:,3) - q1));   % 1사분위수에 가장 가까운 점
[~,t] = min(abs(P_use(:,3) - q3));


% 3) s–t 최소컷 (max-flow)
[~, ~, cut] = maxflow(G, s, t);          % cut: 논리 벡터 (S 집합)

% --- 여기부터 추가 ---
N = size(P_use,1);
cutmask = false(N,1);
cutmask(cut) = true;              % 인덱스 → 논리 마스크로 변환
% ---------------------

% 라벨 벡터 생성 (K=2)
labels = zeros(N,1,'uint32');
labels(cutmask)  = 1;             % S쪽
labels(~cutmask) = 2;                        % T쪽

% 6) 셀 배열로 추출 (clusters{1}, clusters{2})
K = 2;
cmap = lines(1*K);
clusters = cell(1,K);
for kcid = 1:K
    clusters{kcid} = P_use(labels==kcid, :);
end
clusters = clusters(~cellfun('isempty',clusters));
%%
% 시각화용 서브샘플
Eidx = (1:numel(rows))';                  % 유향 간선 인덱스
max_edges_to_draw = 8000;                % 필요시 조정
if numel(Eidx) > max_edges_to_draw
    Eidx = Eidx(randperm(numel(Eidx), max_edges_to_draw));
end

ru = rows(Eidx); cu = cols(Eidx);

X = [P_use(ru,1) P_use(cu,1) NaN(size(ru))];
Y = [P_use(ru,2) P_use(cu,2) NaN(size(ru))];
Z = [P_use(ru,3) P_use(cu,3) NaN(size(ru))];

cmap = parula(256);  % 또는 jet(256)
colormap(cmap);
C = interp1(linspace(0, 1, size(cmap,1)), cmap, w_sim(Eidx), 'linear', 'extrap');

figure(3);  axis equal; box on;
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), 8, [0.8 0.8 0.8], 'filled');hold on;
%line(X', Y', Z', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.4]);
for i = 1: numel(ru)
    line(X(i,:), Y(i,:), Z(i,:), 'Color', C(i,:), 'LineWidth', 0.6);
end
colorbar; caxis([0 1]);
view([1,1,1]); axis equal;;
title(sprintf('kNN links (drawn %d / %d)', numel(ru), numel(rows)));
%%

figure(12); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; 
K = max(labels);
if K > 0
    cmap = lines(K);
    for kcid = 1:K
        idx = (labels == kcid);
        if any(idx)
            scatter3(P_use(idx,1), P_use(idx,2), P_use(idx,3), ...
                     3, cmap(kcid,:), 'filled', ...
                     'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6);
        end
    end
end
hold off
axis equal
xlim([-30.5 30.5]); ylim([-25.5 25.5]); zlim([-35.5 30.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');
%% gomory-hu
T = gomory_hu_tree(G);   % G는 위에서 만든 대칭화+max simplify 그래프

% 트리를 원 좌표 위에 표시(간선 굵기는 컷값에 비례)
[sT,tT,wT] = find(adjacency(T,'weighted'));
lw = 0.2 + 1.0 * (wT - min(wT)) / max(eps, (max(wT)-min(wT)));

X = [P_use(sT,1) P_use(tT,1) NaN(size(sT))];
Y = [P_use(sT,2) P_use(tT,2) NaN(size(sT))];
Z = [P_use(sT,3) P_use(tT,3) NaN(size(sT))];

figure(14); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), 8, 'k', 'filled');
hold on; axis equal; box on;
hl = line(X', Y', Z', 'Color', [0.85 0.2 0.2 0.3]);
hl = hl(:);
for e=1:numel(hl), set(hl(e),'LineWidth',lw(e)); end
view(3); camproj('perspective'); camlight; lighting gouraud;
title('Gomory–Hu tree over point cloud');
hold off
figure(15); 
scatter(P_use(:,1), P_use(:,3), 8, 'k', 'filled'); 
hold on; axis equal; box on;
line(X', Z', 'Color', [0.85 0.2 0.2 0.1]);          % 투명한 선
title('GH tree projected to XY-plane');
hold off