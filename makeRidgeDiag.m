function d = makeRidgeDiag(order, lambda, s)
% MAKERIDGEDIAG  Sobolev 기반 차수별 ridge 페널티 벡터 생성
%
% ── Sobolev 모드 (lambda 가 스칼라) ─────────────────────────────────
%   d = makeRidgeDiag(order, lambda)      % s=1  (H^1 Sobolev)
%   d = makeRidgeDiag(order, lambda, s)   % H^s Sobolev
%
%
%     level l 의 페널티 = lambda × [l(l+1)]^s
%       l=0 : 0×1     =  0         (상수항 → 페널티 없음)
%       l=2 : 2×3     =  6  → ×λ
%       l=4 : 4×5     = 20  → ×λ
%       l=6 : 6×7     = 42  → ×λ
%
%
% ── 공통 출력 ───────────────────────────────────────────────────────
%   d : nT×1 벡터 — regressionFourthOrder 의 lambda 인수로 직접 전달
%       내부에서  A = Phi'Phi + N·diag(d)  로 적용됨
%
% ── 사용 예 ─────────────────────────────────────────────────────────
%   % Sobolev H^1, 전역 lambda=1
%   d = makeRidgeDiag(6, 1);
%   % → l=6: 42, l=4: 20, l=2: 6, l=0: 0
%
%   % Sobolev H^2 (고차 더 강하게)
%   d = makeRidgeDiag(6, 1, 2);
%   % → l=6: 42²=1764, l=4: 400, l=2: 36, l=0: 0
%
%   % 수동 지정
%   d = makeRidgeDiag(6, [0.1, 0.02, 0.003, 0]);
%
%   ransacPar.lambda = d;   % PoliNavigationSolver3_3 에 그대로 전달

    levels = order:-2:0;
    nLev   = numel(levels);

    if nargin < 3 || isempty(s)
        s = 1;
    end

    if isscalar(lambda)
        % ── Sobolev 모드 ──────────────────────────────────────────
        % eigenvalue = l*(l+1), s차 Sobolev → [l*(l+1)]^s
        lambda_per_level = lambda * (levels .* (levels + 1)) .^ s / (order*(order+1))^s;
        % l=0: 0^s = 0 (상수항 자동 제외)

        % 참고 출력
        fprintf('makeRidgeDiag: Sobolev H^%g  (order=%d, λ=%.3g)\n', s, order, lambda);
        for i = 1:nLev
            l = levels(i);
            fprintf('  l=%d : l(l+1)=%2d  →  penalty = %.4g\n', ...
                l, l*(l+1), lambda_per_level(i));
        end
    else
        % ── 수동 모드 ─────────────────────────────────────────────
        if numel(lambda) ~= nLev
            error('makeRidgeDiag: lambda 벡터 길이(%d) != level 수(%d).\n  order=%d → levels=%s', ...
                  numel(lambda), nLev, order, mat2str(levels));
        end
        lambda_per_level = lambda(:).';
    end

    % 항 배치: level l → 2l+1 개 항 (homogeneFischerTerms 와 동일 순서)
    d = [];
    for i = 1:nLev
        l      = levels(i);
        nTerms = 2*l + 1;
        d      = [d; lambda_per_level(i) * ones(nTerms, 1)]; %#ok<AGROW>
    end
end
