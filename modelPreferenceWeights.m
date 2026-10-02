function [Wp, info] = modelPreferenceWeights(P, W, opts)
% 다항식 모델 선호(preference) 커널로 그래프 가중치 재조정
%
% 아이디어: "같은 다항식 모델이 두 점을 함께 설명하면 강하게 잇는다".
% 국소 패치마다 음함수 모델 f_m(x) = sum_k beta_k * phi_k(x - c_m) 을 하나씩 피팅해
% 가설 M개를 만들고, 점마다 가설별 적합도 벡터(선호 벡터)를 만든 뒤
% 두 점의 선호 벡터 유사도(Tanimoto)를 간선 가중치에 곱한다.
% 잔차 기반이라 점 밀도가 불균일해도 영향을 덜 받는다 (T-linkage / RPA 계열).
%
% 입력
%   P    : N x 3 점군 (W와 같은 순서)
%   W    : N x N sparse 유사도 행렬 (buildGraph 출력)
%   opts : .order     다항식 차수 (4 또는 6)            (기본 4)
%          .M         가설 개수                          (기본 200)
%          .kPatch    가설 하나를 피팅할 패치 점 개수     (기본 250)
%          .sigmaD    Taubin 거리 스케일                 (기본 0.3)
%          .lambda    피팅 ridge 계수 (Sobolev)          (기본 1e-3)
%          .beta      Tanimoto 지수 (클수록 강하게 반영)  (기본 1)
%          .seed      난수 시드                          (기본 11)
%
% 출력
%   Wp   : 가중치가 조정된 N x N sparse 행렬 (간선 집합은 W와 동일)
%   info : .dT      N x M Taubin 거리 행렬
%          .nHyp    실제로 성공한 가설 수
%          .T       간선별 Tanimoto 유사도
%          .wi,.wj  대응하는 간선 목록 (상삼각)
%
% 주의: homogeneFischerTerms 가 Symbolic Math Toolbox 를 쓴다.

if nargin < 3, opts = struct(); end
order  = getfield_def(opts, 'order', 4);
M      = getfield_def(opts, 'M', 200);
kPatch = getfield_def(opts, 'kPatch', 250);
sigmaD = getfield_def(opts, 'sigmaD', 0.3);
lambda = getfield_def(opts, 'lambda', 1e-3);
betaEx = getfield_def(opts, 'beta', 1);
seed   = getfield_def(opts, 'seed', 11);

N = size(P,1);
rng(seed);

% 1) 가설 생성: 무작위 시드 주변 패치마다 음함수 다항식 하나
TermsC = homogeneFischerTerms(order);
[Funcs, Grads, nT] = makeFuncsGradsStack(TermsC);
idxP  = knnsearch(P, P, 'K', min(kPatch, N));
seeds = randperm(N, min(M, N));

dT = zeros(N, numel(seeds)); nHyp = 0;
for m = 1:numel(seeds)
    patch = P(idxP(seeds(m),:), :);
    cm    = mean(patch, 1);
    try
        beta = regressionFourthOrder(patch - cm, Funcs, lambda, order, 1);
    catch
        continue
    end
    if any(~isfinite(beta)), continue; end
    Xc = P - cm;
    f  = Funcs(Xc) * beta - 1;
    g  = reshape(Grads(Xc) * beta, N, 3);
    nHyp = nHyp + 1;
    dT(:, nHyp) = abs(f) ./ max(vecnorm(g, 2, 2), eps);   % Taubin 1차 거리
end
dT = dT(:, 1:nHyp);

% 2) 선호 벡터 → 간선별 Tanimoto 유사도
Pref = exp(-(dT / sigmaD).^2);
[wi, wj, wval] = find(triu(W, 1));
ip = sum(Pref(wi,:) .* Pref(wj,:), 2);
n1 = sum(Pref(wi,:).^2, 2);
n2 = sum(Pref(wj,:).^2, 2);
T  = ip ./ max(n1 + n2 - ip, eps);

w  = wval .* (T .^ betaEx);
Wp = sparse([wi; wj], [wj; wi], [w; w], N, N);

info = struct('dT', dT, 'nHyp', nHyp, 'nT', nT, 'T', T, 'wi', wi, 'wj', wj);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
