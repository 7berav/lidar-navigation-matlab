function [labels, info] = cutRecursive(W, opts)
% 원 그래프 위에서 재귀 이분할. 기본은 "안 자름" 이고 병목이 확실할 때만 자른다
% (임베딩 MST 를 만들지 않는다. K 를 미리 정하지 않는다.)
%
% 한 조각에 대해
%   1) 그 조각만의 라플라시안에서 고유벡터 y_2..y_(nVec+1) 을 구한다
%   2) 각 고유벡터 값 순서로 훑으며(sweep) 모든 문턱 t 에서
%          I(t) = cut(t) / ( dbar * sqrt(n_small(t)) )
%      를 계산한다. 표면(2차원) 위에서 "경계 길이 / sqrt(넓이)" 에 해당하는 등주비라
%      점 개수에 무관하다 (conductance 는 점이 많을수록 아무 컷이나 작아진다).
%   3) 병목 정도  rho = min_t I(t) / median_t I(t)
%      균일한 판·띠·구는 길쭉함과 무관하게 rho 가 0.7~1 이고, 목이나 오목 접합은 rho << 1.
%   4) rho <= tau 일 때만 자르고 양쪽을 다시 큐에 넣는다. 아니면 그 조각은 확정.
%
% 정지 규칙은 두 가지 (opts.stop)
%   'rho'   : 위의 등주비 병목 판정 (기본)
%   'ratio' : 부분 그래프의 lambda_3/lambda_2 > r 일 때만 분할 (두 덩어리 구조가 뚜렷할 때).
%             분할 위치는 Fiedler 벡터 sweep 에서 conductance 가 최소인 문턱.
%
% 입력
%   W    : N x N sparse 대칭 유사도 (buildGraphConvex / buildGraph 출력)
%   opts : .stop      'rho' | 'ratio'                           (기본 'rho')
%          .tau       rho 문턱, 작을수록 덜 자른다              (기본 0.35)
%          .r         ratio 문턱, 클수록 덜 자른다              (기본 5)
%          .minSize   조각 최소 점 수 (양쪽 모두)               (기본 max(60, 1% N))
%          .nVec      훑어볼 고유벡터 수                        (기본 3)
%          .alpha     Coifman–Lafon alpha                       (기본 1)
%          .maxPieces 조각 수 상한                              (기본 Inf)
%          .gt        (진단 전용) N x 1 정답 라벨. 분할 결정에는 쓰지 않고
%                     info.nodes(i).gtTrue 에 "이 분할이 부품 경계였는가" 만 기록
%
% 출력
%   labels : N x 1 uint32 (1부터, 0 없음 — 작은 부스러기는 가장 강하게 이어진 조각에 붙임)
%   info   : .nodes  시도한 분할마다 [n, rho, Imin, accepted, gtTrue, vec, lamRatio]

if nargin < 2, opts = struct(); end
N       = size(W,1);
stopRule = getfield_def(opts, 'stop', 'rho');
rThr    = getfield_def(opts, 'r', 5);
tau     = getfield_def(opts, 'tau', 0.35);
minSize = getfield_def(opts, 'minSize', max(60, round(0.01*N)));
nVec    = getfield_def(opts, 'nVec', 3);
alpha   = getfield_def(opts, 'alpha', 1);
maxP    = getfield_def(opts, 'maxPieces', Inf);
gt      = getfield_def(opts, 'gt', []);

labels = zeros(N, 1, 'uint32'); nLab = 0;
nodes  = zeros(0, 7);

comp  = conncomp(graph(W, 'upper')).';
queue = {};
for c = 1:max(comp)
    ii = find(comp == c);
    if numel(ii) >= minSize, queue{end+1} = ii; end %#ok<AGROW>
end

while ~isempty(queue)
    idx = queue{1}; queue(1) = [];                  % 너비 우선 (큰 구조부터)
    n = numel(idx);
    done = true;
    if n >= 2*minSize && (nLab + numel(queue) + 1) < maxP
        [maskA, s] = trySplit(W(idx,idx), minSize, nVec, alpha, stopRule);
        if ~isempty(maskA)
            if strcmpi(stopRule, 'ratio'), acc = s.lamRatio > rThr; else, acc = s.rho <= tau; end
            nodes(end+1,:) = [n, s.rho, s.Imin, acc, gtCheck(gt, idx, maskA, minSize), s.vec, s.lamRatio]; %#ok<AGROW>
            if acc
                done = false;
                for side = {idx(maskA), idx(~maskA)}
                    jj = side{1};
                    cc = conncomp(graph(W(jj,jj), 'upper')).';
                    for c = 1:max(cc)
                        kk = jj(cc == c);
                        if numel(kk) >= minSize, queue{end+1} = kk; end %#ok<AGROW>
                    end
                end
            end
        end
    end
    if done
        nLab = nLab + 1; labels(idx) = nLab;
    end
end

labels = attachLeftover(labels, W);
info = struct('nodes', array2table(nodes, 'VariableNames', ...
              {'n','rho','Imin','accepted','gtTrue','vec','lamRatio'}), 'minSize', minSize, 'tau', tau);
end

function [maskA, s] = trySplit(Ws, minSize, nVec, alpha, stopRule)
% 연결된 조각 하나에서 가장 병목다운 문턱 분할을 찾는다
maskA = []; s = struct('rho', NaN, 'Imin', NaN, 'vec', NaN, 'lamRatio', NaN);
n = size(Ws,1);
d = full(sum(Ws,2)); dbar = mean(d);

[L, Dis] = normalizeGraph(Ws, alpha);
kE = min(max(nVec + 1, 3), n - 2);
try
    [U, lam] = eigs(L + 1e-9*speye(n), kE, 'smallestabs', ...
                    'Tolerance', 1e-5, 'MaxIterations', 500, 'Display', 0);
catch
    return
end
[lam, o] = sort(diag(lam) - 1e-9, 'ascend'); U = U(:,o);
Y = Dis * U(:, 2:min(nVec+1, kE));
lamRatio = NaN;
if kE >= 3, lamRatio = lam(3) / max(lam(2), 1e-12); end

[ei, ej, ew] = find(triu(Ws, 1));
t  = (1:n-1)';
nS = min(t, n - t);
valid = nS >= minSize;
if ~any(valid), return; end

best = Inf;
for v = 1:size(Y,2)
    [ysort, ord] = sort(Y(:,v));
    pos = zeros(n,1); pos(ord) = 1:n;
    assocPre = cumsum(accumarray(max(pos(ei), pos(ej)), ew, [n 1]));
    volPre   = cumsum(d(ord));
    cutPre   = volPre(1:n-1) - 2*assocPre(1:n-1);        % 앞쪽 t 개와 나머지 사이의 컷
    I = cutPre ./ (dbar * sqrt(nS));
    I(~valid) = Inf;
    [Imin, tStar] = min(I);
    rho = Imin / median(I(valid));
    if strcmpi(stopRule, 'ratio')                        % Fiedler 벡터의 최소 conductance 문턱
        phi = cutPre ./ min(volPre(1:n-1), volPre(n) - volPre(1:n-1));
        phi(~valid) = Inf;
        [~, tStar] = min(phi);
    end
    if rho < best
        best  = rho;
        maskA = false(n,1); maskA(ord(1:tStar)) = true;
        s = struct('rho', rho, 'Imin', Imin, 'vec', v + 1, 'lamRatio', lamRatio, 'thr', ysort(tStar));
    end
    if strcmpi(stopRule, 'ratio'), break; end            % ratio 규칙은 Fiedler 벡터만 사용
end
end

function ok = gtCheck(gt, idx, maskA, minSize)
% (진단) 분할이 부품 경계와 일치하는가: 주요 부품이 모두 한쪽에 90% 이상 몰리고,
% 양쪽에 서로 다른 주요 부품이 하나 이상 있다. gt 가 없으면 NaN
ok = NaN;
if isempty(gt), return; end
g = double(gt(idx));
labs = unique(g);
cntA = arrayfun(@(l) nnz(g == l &  maskA), labs);
cntB = arrayfun(@(l) nnz(g == l & ~maskA), labs);
major = (cntA + cntB) >= minSize/2;
fa = cntA(major) ./ (cntA(major) + cntB(major));
ok = double(all(fa <= 0.1 | fa >= 0.9) && any(fa >= 0.9) && any(fa <= 0.1));
end

function lab = attachLeftover(lab, W)
% 라벨 0 (minSize 미만 부스러기)을 W 가중치 합이 가장 큰 이웃 조각에 붙임
lab = double(lab); K = max(lab);
if K == 0, lab(:) = 1; lab = uint32(lab); return; end
for it = 1:100
    un = find(lab == 0);
    if isempty(un), break; end
    [ii, jj, ww] = find(W(un, :));
    ii = ii(:); jj = jj(:); ww = ww(:);           % un 이 1개면 find 가 행벡터를 돌려준다
    ok = lab(jj) > 0;
    if ~any(ok), break; end
    A = accumarray([ii(ok), lab(jj(ok))], ww(ok), [numel(un) K]);
    [mx, bestK] = max(A, [], 2);
    upd = mx > 0;
    if ~any(upd), break; end
    lab(un(upd)) = bestK(upd);
end
if any(lab == 0)                                   % 어디에도 안 이어진 고립 조각
    lab(lab == 0) = K + 1;
end
lab = uint32(lab);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
