function [K, info] = selectKRatio(lambda, opts)
% 고유값 "비율" 간격으로 군집 수 K 선택 (selectEigenGap 의 차이 간격과 달리 스케일 불변)
%
%   K1  'max'      : K = argmax_{K>=2} lambda_{K+1} / lambda_K
%   K2  'first'    : 비율이 r0 를 넘는 가장 작은 K
%   K1s 'skipnear' : "거의 0" 인 고유값 nNear 개(끊어졌거나 가는 연결로 거의 끊어진 조각)는
%                    각각 별도 군집으로 확정하고, 그 뒤 구조적 고유값에서만 최대 비율을 찾는다.
%                    K = argmax_{k > nNear} lambda_{k+1}/lambda_k,  단 그 비율이 rMin 미만이면
%                    더 나눌 구조가 없다고 보고 K = nNear (nNear = 1 이면 K = 1 도 가능).
%                    "거의 0" = lambda < nearRel * lambda_(kMax+1)  (스케일에 맞춰 상대적으로)
%
% 'max'/'first' 에서는 floor 보다 작은 고유값을 수치적으로 0 (= 연결성분) 으로 본다. 영고유값이
% n0 >= 2 개면 그 뒤의 간격이 사실상 무한대이므로 K = n0 를 돌려준다 (끊어진 덩어리 = 군집).
%
% 입력
%   lambda : 고유값 오름차순 (embedSpectral 출력)
%   opts   : .rule    'max' | 'first' | 'skipnear'   (기본 'max')
%            .r0      'first' 의 문턱                 (기본 2)
%            .floor   수치 잡음 하한                  (기본 1e-7)
%            .kMax    탐색 상한                       (기본 min(15, numel(lambda)-1))
%            .nearRel 'skipnear': 거의 0 판정 비율    (기본 0.01)
%            .rMin    'skipnear': 추가 분할에 필요한 최소 비율 (기본 2)
%
% 출력
%   K    : 추정 군집 수
%   info : .ratio (ratio(k) = lambda_{k+1}/lambda_k), .n0 (영고유값 수), .nNear, .ratioAtK

if nargin < 2, opts = struct(); end
rule = getfield_def(opts, 'rule', 'max');
r0   = getfield_def(opts, 'r0', 2);
fl   = getfield_def(opts, 'floor', 1e-7);
kMax = getfield_def(opts, 'kMax', min(15, numel(lambda)-1));

lam   = max(lambda(:), fl);
ratio = lam(2:end) ./ lam(1:end-1);
n0    = max(nnz(lambda(:) < fl), 1);
kMax  = min(kMax, numel(ratio));
nNear = NaN;

if strcmpi(rule, 'skipnear')
    nearRel = getfield_def(opts, 'nearRel', 0.01);
    rMin    = getfield_def(opts, 'rMin', 2);
    lamRef  = lambda(min(kMax+1, numel(lambda)));
    nNear   = max(nnz(lambda(1:kMax) < nearRel * lamRef), 1);
    ks = (nNear+1):kMax;
    K  = nNear;
    if ~isempty(ks)
        [rm, i] = max(ratio(ks));
        if rm >= rMin, K = ks(i); end
    end
elseif n0 >= 2
    K = n0;
else
    ks = 2:kMax;
    switch lower(rule)
        case 'max'
            [~, i] = max(ratio(ks)); K = ks(i);
        case 'first'
            i = find(ratio(ks) > r0, 1);
            if isempty(i), [~, i] = max(ratio(ks)); end      % 넘는 곳이 없으면 최대 비율
            K = ks(i);
        otherwise
            error('selectKRatio: 알 수 없는 rule "%s"', rule);
    end
end
info = struct('ratio', ratio, 'n0', n0, 'nNear', nNear, 'ratioAtK', ratio(max(K,1)));
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
