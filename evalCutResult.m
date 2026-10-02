function r = evalCutResult(lab, gt, Wt, lam)
% 분할 결과 채점 (정답 기준 + 그래프 기준)
%
% 입력
%   lab : N x 1 예측 라벨 (0 은 별도 군집으로 취급)
%   gt  : N x 1 정답 부품 라벨
%   Wt  : N x N 그래프 (NCut·phi 계산용, 대칭)
%   lam : 고유값 오름차순 (Upsilon 계산용, 생략 가능)
%
% 출력 r (struct)
%   .K .ARI .purity(부품을 덜 합칠수록 1) .frag(부품의 90% 를 덮는 군집 수 평균, 1 이 정답)
%   .ncutK (K-way NCut)  .upsilon (lambda_{K+1}/max phi)  .minVolFrac  .maxPhi  .labels (1..K)

if nargin < 4, lam = []; end
lab = double(lab(:));
if any(lab == 0), lab(lab == 0) = max(lab) + 1; end
[~, ~, li] = unique(lab); [gu, ~, gi] = unique(double(gt(:)));
K  = max(li);
e  = evalSegmentation(uint32(li), gt);
Ct = accumarray([li gi], 1);
fr = zeros(numel(gu), 1);
for q = 1:numel(gu)
    col = sort(Ct(:,q), 'descend'); fr(q) = find(cumsum(col) >= 0.9*sum(col), 1);
end
major = sum(Ct, 1).' >= 30;

d = full(sum(Wt,2)); volTot = sum(d);
vol = accumarray(li, d, [K 1]);
[ei, ej, ew] = find(Wt);
in = li(ei) == li(ej);
cutv = vol - accumarray(li(ei(in)), ew(in), [K 1]);
phi  = cutv ./ max(min(vol, volTot - vol), eps);
ups  = NaN;
if K >= 2 && K + 1 <= numel(lam), ups = lam(K+1) / max(max(phi), eps); end

r = struct('K', K, 'ARI', e.ARI, 'purity', sum(max(Ct,[],2)) / sum(Ct(:)), 'frag', mean(fr(major)), ...
           'ncutK', sum(cutv ./ max(vol, eps)), 'upsilon', ups, 'minVolFrac', min(vol)/volTot, ...
           'maxPhi', max(phi), 'labels', uint32(li));
end
