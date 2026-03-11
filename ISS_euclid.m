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

% 2) 반경 τ 이웃 리스트(셀 배열) — 범위 탐색
tau = 0.40;  % [좌표 단위] 데이터 스케일에 맞게
nbrs = rangesearch(Mdl, P_use, tau);
for i = 1:numel(nbrs)
    nbrs{i}(nbrs{i}==i) = [];   % self index 제거
end
% 3) 연결 성분 라벨링 (BFS/Union-Find; 그래프 굳이 만들지 않고 한 번만 방문)
N = size(P_use,1);
labels = zeros(N,1,'uint32');
cid = uint32(0);
visited = false(N,1);

%%
for i = 1:N
    if visited(i), continue; end
    cid = cid + 1;
    % BFS
    q = i;
    visited(i) = true;
    labels(i) = cid;
    head = 1;
    while head <= numel(q)
        u = q(head); head = head + 1;
        nu = nbrs{u};
        % 자기 자신 제외
        nu = nu(nu ~= u);
        % 아직 방문 안 한 점들만
        nu = nu(~visited(nu));
        visited(nu) = true;
        labels(nu) = cid;
        % 큐 확장
        q = [q; nu(:)]; %#ok<AGROW>
    end
end
%%
% 4) 너무 작은 조각 제거
numel(unique(labels))
Nmin = 60;  % RANSAC 최소 표본보다 충분히 크게
[uc,~,ic] = unique(labels);
counts = accumarray(ic,1);
keep = ismember(labels, uc(counts>=Nmin));
labels(~keep) = 0;  % 0 = 버림
numel(unique(labels))
% 5) 셀 배열로 클러스터 추출
K = double(max(labels));
clusters = cell(1,K);
for k = 1:K
    clusters{k} = P_use(labels==k,:);
end
clusters = clusters(~cellfun('isempty',clusters));

%%
figure(12); 
scatter3(P_use(:,1), P_use(:,2), P_use(:,3), ...
         1, [0.5 0.5 0.5], 'filled', ...
         'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
hold on; axis equal
cmap = lines(numel(clusters));
for i = 1:numel(clusters)
    Ci = clusters{i};
    scatter3(Ci(:,1),Ci(:,2),Ci(:,3), 3, cmap(i,:), 'filled', ...
                     'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.9);
end
hold off
axis equal
xlim([-30.5 30.5]); ylim([-10.5 15.5]); zlim([-25.5 20.5]);   % 필요에 맞게
xlabel('X'); ylabel('Y'); zlabel('Z');
