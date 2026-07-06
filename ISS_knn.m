% 분할 
% kd tree -> knn

P = readmatrix('ISS_stationary.xyz', 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21,:);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8,:);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 &  P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13,:);

P_use = [P1 ; P2; P3; P4];
% 1) KD-tree 생성 (Statistics/ML Toolbox)
Mdl = createns(P_use, 'NSMethod','kdtree', 'Distance','euclidean');

k = 10;                                % 최근접 이웃 수
[idx, dist] = knnsearch(Mdl, P_use, 'K', k+1);
idx  = idx(:,2:end);                    % 자기 자신 제거
dist = dist(:,2:end);

% kNN 간선 리스트 (유향) → 그래프 만들 때 대칭화 예정
rows = repmat((1:size(P_use,1))', k,1);
cols = idx(:);
w    = dist(:);  
%%
% 1) 무향 그래프(가중: 거리)
G = graph(rows, cols, w, size(P_use,1));
G = simplify(G, 'min');                 % 중복 간선 있으면 최소 가중치만

% 2) 최소 신장 트리(MST)
T = minspantree(G, 'Method','sparse');

% === 옵션 1: 가장 긴 간선 C개 제거 ===
 C = 5;                                  % 자를 간선 개수 (작게 시작)
[~, ord] = sort(T.Edges.Weight, 'descend');
cut_idx = ord(1:min(C, numedges(T)));
T2 = rmedge(T, T.Edges.EndNodes(cut_idx,1), T.Edges.EndNodes(cut_idx,2));

% % === 옵션 2: 길이 임계 이상 제거 ===
%T_cut = 0.3;                          % 예) 컷 임계값(데이터 스케일 맞춰 조정)
%long_e = find(T.Edges.Weight >= T_cut);
%T2 = rmedge(T, T.Edges.EndNodes(long_e,1), T.Edges.EndNodes(long_e,2));

% 3) 연결 성분 라벨
labels = conncomp(T2)';                 % 1..K

% 4) 너무 작은 성분 제거(옵션)
Nmin = 40;
K = max(labels);
numel(unique(labels))
counts = accumarray(labels(:),1,[K 1]);
small = ismember(labels, find(counts < Nmin));
labels(small) = 0;                      % 0 = 버림
numel(unique(labels))

% 5) (선택) 셀 배열로 추출
K = max(labels);
clusters = cell(1,K);
for kcid = 1:K
    clusters{kcid} = P_use(labels==kcid, :);
end
clusters = clusters(~cellfun('isempty',clusters));

%%
figure(12); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; 
K = max(labels);
if K > 0
    cmap = lines(K+10);
    for kcid = 1:K
        idx = (labels == kcid);
        if any(idx)
            scatter3(P_use(idx,1), P_use(idx,2), P_use(idx,3), ...
                     2, cmap(kcid+1,:), 'filled', ...
                     'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6);
        end
    end
end
hold off
axis equal
xlim([-30.5 30.5]); ylim([-15.5 15.5]); zlim([-25.5 20.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');
