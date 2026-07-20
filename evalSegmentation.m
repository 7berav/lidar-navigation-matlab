function metrics = evalSegmentation(labels, labels_gt)
% 분할 결과 정량 평가 — 라벨 순열 불변 지표 (ARI)
%
% 입력
%   labels    : N x 1 예측 라벨 (0 = 버려진 점, 평가에서 제외)
%   labels_gt : N x 1 정답 라벨
%
% 출력 metrics:
%   .ARI      Adjusted Rand Index (1=완벽, 0=무작위 수준)
%   .coverage 라벨이 부여된 점 비율 (labels>0)
%   .K_pred   예측 군집 수 / .K_gt 정답 군집 수
%
% 주의: 라벨 번호 직접 비교 금지 (같은 분할도 순열이 다르면 틀렸다고 나옴)

labels    = double(labels(:));
labels_gt = double(labels_gt(:));
assert(numel(labels) == numel(labels_gt), 'evalSegmentation: 길이 불일치');

valid = labels > 0;
metrics.coverage = mean(valid);

la = labels(valid);
lb = labels_gt(valid);

% contingency table
[~,~,ia] = unique(la);
[~,~,ib] = unique(lb);
Ct = accumarray([ia, ib], 1);

nij = sum(sum(Ct .* (Ct-1) / 2));
ai  = sum(Ct, 2); a = sum(ai .* (ai-1) / 2);
bj  = sum(Ct, 1); b = sum(bj .* (bj-1) / 2);
n   = numel(la);
nC2 = n*(n-1)/2;

expected = a * b / nC2;
maxIdx   = (a + b) / 2;
if maxIdx - expected < eps
    metrics.ARI = 1;    % 퇴화 케이스 (양쪽 모두 군집 1개 등)
else
    metrics.ARI = (nij - expected) / (maxIdx - expected);
end

metrics.K_pred = numel(unique(la));
metrics.K_gt   = numel(unique(lb));
end
