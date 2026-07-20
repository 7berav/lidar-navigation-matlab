function [labels, diag_out] = segmentEmbedding(Yembed, W, opts)
% 고유벡터 임베딩 공간에서 MST 기반 분할
% (ISS_normcut_4_modified.m 294~426행의 함수화 + conductance 판정 추가)
%
% 입력
%   Yembed : N x m 임베딩 좌표 (예: DinvSqrt*U 의 열 일부)
%   W      : N x N 원 그래프 유사도 행렬 (conductance 판정용)
%   opts   : .ky        임베딩 kNN 이웃 수               (기본 40)
%            .cutMethod 'length' | 'conductance'         (기본 'conductance')
%            .qThr      length: 컷 분위수 / conductance: 후보 분위수 (기본 0.96)
%            .phiMax    conductance 채택 임계값          (기본 0.05)
%            .minSize   최소 군집 크기, 미달 라벨=0      (기본 60)
%
% 출력
%   labels   : N x 1 uint32 (0 = 버려진 점)
%   diag_out : .cutEdges  실제 잘린 간선 [u v len phi]
%              .candEdges 후보 간선 [u v len]
%              .T_all     최종 포레스트 graph 객체
%
% cutMethod 이론 근거 (가이드 §4.5):
%   'length'      : MST 긴 간선 절단 (현행 baseline, 사슬에 취약)
%   'conductance' : 후보 간선 제거 시 원 그래프 W에서의 conductance
%                   phi = cut(S,S~)/min(vol S, vol S~) 를 계산해
%                   phi <= phiMax 인 컷만 채택 (Cheeger 근거화)

if nargin < 3, opts = struct(); end
ky      = getfield_def(opts, 'ky', 40);
method  = getfield_def(opts, 'cutMethod', 'conductance');
qThr    = getfield_def(opts, 'qThr', 0.96);
phiMax  = getfield_def(opts, 'phiMax', 0.05);
minSize = getfield_def(opts, 'minSize', 60);

N = size(Yembed,1);

% --- 1) 임베딩 공간 kNN 그래프 → 성분별 MST (포레스트) ---
Mdl = createns(Yembed, 'NSMethod','kdtree', 'Distance','euclidean');
[idxNN, distNN] = knnsearch(Mdl, Yembed, 'K', ky+1);
idxNN  = idxNN(:,2:end);
distNN = distNN(:,2:end);

I = repmat((1:N)', ky, 1);
J = idxNN(:);
V = distNN(:);
S = sparse(I, J, V, N, N);
S = max(S, S.');
G_knn = graph(S, 'upper');

comp = conncomp(G_knn);
m = max(comp);

% 포레스트 간선을 전역 인덱스로 수집
E_global = zeros(0,2);
W_global = zeros(0,1);
for i = 1:m
    Vi = find(comp == i);
    Ti = minspantree(subgraph(G_knn, Vi), 'Method','sparse');
    E  = Ti.Edges.EndNodes;
    E_global = [E_global; [Vi(E(:,1))', Vi(E(:,2))']]; %#ok<AGROW>
    W_global = [W_global; Ti.Edges.Weight];            %#ok<AGROW>
end

% --- 2) 컷 후보: 분위수 이상 긴 간선, 긴 것부터 ---
tauLen = quantile(W_global, qThr);
candMask = W_global >= tauLen;
candIdx  = find(candMask);
[~, ord] = sort(W_global(candIdx), 'descend');
candIdx  = candIdx(ord);
candEdges = [E_global(candIdx,:), W_global(candIdx)];

% 포레스트 인접 리스트 (컷 반영하며 갱신)
adjF = sparse(E_global(:,1), E_global(:,2), 1:size(E_global,1), N, N);
adjF = adjF + adjF.';
removed = false(size(E_global,1),1);

dW = full(sum(W,2));            % 원 그래프 차수 (conductance 볼륨용)
volTotal = sum(dW);

cutEdges = zeros(0,4);
for c = 1:numel(candIdx)
    e = candIdx(c);
    u = E_global(e,1); v = E_global(e,2);

    switch lower(method)
        case 'length'
            phi = NaN;
            accept = true;
        case 'conductance'
            % 간선 e를 제거했을 때 u쪽 서브트리 S를 BFS로 수집
            Smask = subtreeMask(adjF, removed, u, e, N);
            volS  = sum(dW(Smask));
            volSc = volTotal - volS;
            % cut(S, S~) = vol(S) - 2*assoc(S,S) ... 대신 직접 계산
            cutVal = sum(sum(W(Smask, ~Smask)));
            phi = cutVal / max(min(volS, volSc), eps);
            accept = (phi <= phiMax);
        otherwise
            error('segmentEmbedding: 알 수 없는 cutMethod "%s"', method);
    end

    if accept
        removed(e) = true;
        cutEdges = [cutEdges; u, v, W_global(e), phi]; %#ok<AGROW>
    end
end

% --- 3) 라벨: 컷 반영된 포레스트의 연결 성분 ---
Ekeep = E_global(~removed, :);
T_all = graph(Ekeep(:,1), Ekeep(:,2), W_global(~removed), N);
labels_all = uint32(conncomp(T_all))';

% --- 4) 최소 크기 필터 → 연속 라벨 재맵핑, 미달=0 ---
K_now = double(max(labels_all));
sz = accumarray(double(labels_all), 1, [K_now, 1]);
keep = find(sz >= minSize);
labels = zeros(N,1,'uint32');
if ~isempty(keep)
    remap = zeros(K_now,1,'uint32'); remap(keep) = 1:numel(keep);
    labels = remap(labels_all);
end

diag_out = struct('cutEdges', cutEdges, 'candEdges', candEdges, ...
                  'T_all', T_all);
end

function Smask = subtreeMask(adjF, removed, src, eSkip, N)
% 포레스트에서 간선 eSkip을 막은 채 src에서 도달 가능한 노드 집합
Smask = false(N,1);
stack = src;
Smask(src) = true;
while ~isempty(stack)
    u = stack(end); stack(end) = [];
    [~, nbrs, eids] = find(adjF(u,:));
    for t = 1:numel(nbrs)
        v = nbrs(t); eid = eids(t);
        if eid == eSkip || removed(eid) || Smask(v), continue; end
        Smask(v) = true;
        stack(end+1) = v; %#ok<AGROW>
    end
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
