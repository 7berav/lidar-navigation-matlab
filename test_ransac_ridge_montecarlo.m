%% test_ransac_ridge_montecarlo.m
% Ridge vs OLS 단일 드로우 성공률 Monte Carlo 분석
%
% 핵심 아이디어:
%   freezeW=true, freezeOmega=true, freezeIter=true 로 설정
%   → omega 균일 고정 → 매 iter 독립 iid 랜덤 드로우
%   → Log.preScore (Sj) 분포 = 단일 드로우 품질 샘플
%   updateThresh=0 → 로컬 최적화 매 iter 실행 → recall 항상 기록 (iid)
%
% 분석 목표:
%   점수 분포 (preScore, recall) 비교
%   상위 비율 P(Sj > α×w) for α ∈ {0.70, 0.80, 0.90}  [w-스케일 기준]
%   implied N_iter = log(1-p)/log(1-P_success) vs α×w threshold

%% ---- 공통 세팅 ----
totalN      = 10000;
eps_list    = [0.05, 0.10, 0.20];
k_list      = [1.4];
order_list  = [6];
Niter       = 10000;
lambda_list = [0, 1e-3, 3e-3, 1e-2, 3e-2];
conf        = 0.95;

% w-스케일 threshold용 alpha 범위
alpha_arr     = linspace(0.40, 1.02, 100);   % thresh = alpha * w_true
alpha_ref     = [0.70, 0.80, 0.90];          % 상위 비율 bar chart용 기준점

% 데이터 풀 생성
PP11 = generateRandomPointsOnHexagonPrism(totalN + 1000) + randn(totalN + 1000, 3) * 0.007;
PP14 = (2 * rand(5000, 3) - 1) * 2;

%% ---- RANSAC 파라미터 ----
ransacPar = struct( ...
    'maxIter',        Niter,  ...
    'conf',           conf,   ...
    'thresh',         0.40,   ...
    'minInlierRatio', 0.65,   ...
    'updateThresh',   0,      ...  % 0 → 매 iter 로컬 opt 실행 → recall 항상 기록
    'locIters',       4,      ...
    'momentum',       0.66,   ...
    'damping',        0.85,   ...
    'reg',            1e-6,   ...
    'k',              1.4,    ...
    'lambda',         0       ...
);
ransacPar.mc = struct( ...
    'on',          true,    ...
    'saveVarName', '',      ...
    'saveMatFile', '',      ...
    'time',        true,    ...
    'localOff',    false,   ...
    'metric',      'basic'  ...
);
% ★ 핵심: 모든 적응 메커니즘 동결 → 매 iter 독립 iid 드로우
ransacPar.sim = struct( ...
    'freezeW',     true,  ...
    'freezeIter',  true,  ...
    'freezeOmega', true   ...
);

warning('off', 'MATLAB:nearlySingularMatrix');

%% ---- 실행 루프 ----
RES = struct([]);
row = 0;

for order = order_list
    TermsC = homogeneFischerTerms(order);
    [Funcs1, Grads1, ~] = makeFuncsGradsStack(TermsC);
    nT = numel(TermsC);

    for e = eps_list
        n_out   = round(totalN * e);
        n_in    = totalN - n_out;
        if n_in > size(PP11, 1) || n_out > size(PP14, 1)
            warning('데이터 풀 부족: order=%d, e=%.2f → 건너뜀', order, e);
            continue;
        end
        PPm_use = [PP11(1:n_in, :); PP14(1:n_out, :)];
        gtMask  = false(totalN, 1);
        gtMask(1:n_in) = true;
        ransacPar.mc.gtMask = gtMask;

        w_true = n_in / totalN;   % ≈ 1 - eps (ground truth)

        for k_mult = k_list
            ransacPar.k = k_mult;
            k_samples   = round(k_mult * nT);

            for lam = lambda_list
                ransacPar.lambda = lam;

                if lam == 0
                    lam_tag = 'OLS';
                else
                    lam_tag = strrep(sprintf('%.0e', lam), '-0', 'm');
                end
                ransacPar.mc.saveVarName = sprintf( ...
                    'MC_ridge_w%02d_k%02d_or%d_%s', ...
                    round(100*e), round(10*k_mult), order, lam_tag);

                fprintf('order=%d  eps=%.2f  k=%.1f(=%d)  lambda=%-8s ... ', ...
                    order, e, k_mult, k_samples, lam_tag);
                t0 = tic;

                [~, ~, ~, Log] = PoliNavigationSolver3_3_FischerRansac_MC( ...
                    0, PPm_use, order, nT, Funcs1, Grads1, ransacPar);

                elapsed = toc(t0);
                fprintf('done (%.1fs, %d iters)\n', elapsed, numel(Log));

                % ---- 로그 벡터 추출 ----
                preScores  = double([Log.preScore]');
                postScores = double([Log.postScore]');
                recalls    = double([Log.recall]');
                inlierRs   = double([Log.inlierR]');
                F1s        = double([Log.F1]');
                wLog       = double([Log.w]');

                % ---- 결과 저장 ----
                row = row + 1;
                RES(row).order      = order;
                RES(row).eps        = e;
                RES(row).k_mult     = k_mult;
                RES(row).k          = k_samples;
                RES(row).nT         = nT;
                RES(row).lambda     = lam;
                RES(row).lam_tag    = lam_tag;
                RES(row).w_true     = w_true;
                RES(row).wk_theory  = w_true ^ k_samples;
                RES(row).N_theory   = ceil(log(1 - conf) / log(max(realmin, 1 - w_true^k_samples)));
                RES(row).preScores  = preScores;
                RES(row).postScores = postScores;
                RES(row).recalls    = recalls;
                RES(row).inlierRs   = inlierRs;
                RES(row).F1s        = F1s;
                RES(row).wLog       = wLog;
            end
        end
    end
end

%% ---- 분석: w-스케일 threshold sweep ----
% thresh = alpha * w_true  →  eps가 달라도 같은 기준으로 비교 가능
for ri = 1:numel(RES)
    w  = RES(ri).w_true;
    ps = RES(ri).preScores;
    rc = RES(ri).recalls;

    thresh_sj  = alpha_arr * w;   % 절대 Sj 기준 (w-스케일)
    thresh_rec = alpha_arr;       % recall은 0~1 절대값 그대로 사용

    valid_ps = ps(isfinite(ps));
    valid_rc = rc(isfinite(rc));

    P_sj = arrayfun(@(t) mean(valid_ps > t), thresh_sj);
    P_rc = arrayfun(@(t) mean(valid_rc > t), thresh_rec);

    % implied N_iter
    N_impl_sj = ceil(log(1 - conf) ./ log(max(realmin, 1 - P_sj)));

    % alpha_ref 기준 상위 비율 (bar chart용)
    P_ref_sj  = arrayfun(@(a) mean(valid_ps > a * w),  alpha_ref);
    P_ref_rec = arrayfun(@(a) mean(valid_rc > a),       alpha_ref);

    RES(ri).alpha_arr  = alpha_arr;
    RES(ri).thresh_sj  = thresh_sj;
    RES(ri).thresh_rec = thresh_rec;
    RES(ri).P_sj       = P_sj;
    RES(ri).P_rc       = P_rc;
    RES(ri).N_impl_sj  = N_impl_sj;
    RES(ri).P_ref_sj   = P_ref_sj;   % (1 x numel(alpha_ref))
    RES(ri).P_ref_rec  = P_ref_rec;
end

%% ---- 시각화 ----
colors = lines(numel(lambda_list));
lambda_labels = cell(1, numel(lambda_list));
for li = 1:numel(lambda_list)
    if lambda_list(li) == 0
        lambda_labels{li} = 'OLS (\lambda=0)';
    else
        lambda_labels{li} = sprintf('\\lambda=%.0e', lambda_list(li));
    end
end

figBase = 20;

for order = order_list
    for e = eps_list(2)
        for k_mult = k_list

            idx_all = find( ...
                [RES.order]  == order & ...
                abs([RES.eps]    - e)      < 1e-9 & ...
                abs([RES.k_mult] - k_mult) < 1e-9);
            if isempty(idx_all), continue; end

            nT_val    = RES(idx_all(1)).nT;
            k_val     = RES(idx_all(1)).k;
            w_true    = RES(idx_all(1)).w_true;
            wk_theory = RES(idx_all(1)).wk_theory;
            N_theory  = RES(idx_all(1)).N_theory;
            tag       = sprintf('or%d_e%02d_k%02d', order, round(100*e), round(10*k_mult));

            % ---- Fig A: 점수 분포 (preScore + recall 2-panel) ----
            figure(figBase); clf; figBase = figBase + 1;
            set(gcf, 'Position', [50 50 720 560]);

            subplot(2, 1, 1); hold on; grid on;
            for li = 1:numel(lambda_list)
                ri_m = idx_all(li);
                ps = RES(ri_m).preScores;
                histogram(ps(isfinite(ps)), 60, ...
                    'FaceAlpha', 0.35, 'EdgeAlpha', 0.15, ...
                    'FaceColor', colors(li, :), 'DisplayName', lambda_labels{li});
            end
            % alpha_ref 수직선
            for ai = 1:numel(alpha_ref)
                xline(alpha_ref(ai) * w_true, 'k:', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%.0f%%×w', alpha_ref(ai)*100));
            end
            xline(w_true, 'r-', 'LineWidth', 1.5, 'DisplayName', 'w (이론 최대)');
            xlabel('preScore (Sj)'); ylabel('count');
            title(sprintf('preScore 분포  |  order=%d  ε=%.2f  k=%.1f×%d=%d', ...
                order, e, k_mult, nT_val, k_val));
            legend('Location', 'northeast', 'FontSize', 7);
            hold off;

            subplot(2, 1, 2); hold on; grid on;
            for li = 1:numel(lambda_list)
                ri_m = idx_all(li);
                rc = RES(ri_m).recalls;
                histogram(rc(isfinite(rc)), 60, ...
                    'FaceAlpha', 0.35, 'EdgeAlpha', 0.15, ...
                    'FaceColor', colors(li, :), 'DisplayName', lambda_labels{li});
            end
            for ai = 1:numel(alpha_ref)
                xline(alpha_ref(ai), 'k:', 'LineWidth', 1.2);
            end
            xlabel('recall'); ylabel('count');
            title(sprintf('recall 분포  |  updateThresh=0 → 전 iter 로컬 opt'));
            legend('Location', 'northwest', 'FontSize', 7);
            hold off;

            exportgraphics(gcf, sprintf('ridge_dist_%s.png', tag), 'Resolution', 200);

            % ---- Fig B: 상위 비율 bar chart (P(Sj > alpha_ref * w)) ----
            figure(figBase); clf; figBase = figBase + 1;
            set(gcf, 'Position', [800 50 720 400]);

            nAlpha = numel(alpha_ref);
            nLam   = numel(lambda_list);
            barData_sj  = zeros(nAlpha, nLam);
            barData_rec = zeros(nAlpha, nLam);
            for li = 1:nLam
                ri_m = idx_all(li);
                barData_sj(:, li)  = RES(ri_m).P_ref_sj';
                barData_rec(:, li) = RES(ri_m).P_ref_rec';
            end

            subplot(1, 2, 1);
            b = bar(barData_sj);
            for li = 1:nLam
                b(li).FaceColor = colors(li, :);
                b(li).DisplayName = lambda_labels{li};
            end
            set(gca, 'XTickLabel', arrayfun(@(a) sprintf('%.0f%%×w', a*100), ...
                alpha_ref, 'UniformOutput', false));
            ylabel('P(Sj > threshold)');
            title('preScore 상위 비율');
            legend('Location', 'northeast', 'FontSize', 7);
            grid on; ylim([0 1]);

            subplot(1, 2, 2);
            b = bar(barData_rec);
            for li = 1:nLam
                b(li).FaceColor = colors(li, :);
                b(li).DisplayName = lambda_labels{li};
            end
            set(gca, 'XTickLabel', arrayfun(@(a) sprintf('recall>%.2f', a), ...
                alpha_ref, 'UniformOutput', false));
            ylabel('P(recall > threshold)');
            title('recall 상위 비율');
            legend('Location', 'northeast', 'FontSize', 7);
            grid on; ylim([0 1]);

            sgtitle(sprintf('상위 비율  |  ε=%.2f  k=%d  nT=%d', e, k_val, nT_val));
            exportgraphics(gcf, sprintf('ridge_toprate_%s.png', tag), 'Resolution', 200);

            % ---- Fig C: implied N_iter vs alpha (N_theory 점선 제거) ----
            figure(figBase); clf; figBase = figBase + 1;
            set(gcf, 'Position', [800 50 700 420]);
            hold on; grid on;
            for li = 1:nLam
                ri_m = idx_all(li);
                semilogy(RES(ri_m).alpha_arr, RES(ri_m).N_impl_sj, '-', ...
                    'Color', colors(li, :), 'LineWidth', 1.8, ...
                    'DisplayName', lambda_labels{li});
            end
            % alpha_ref 수직선
            for ai = 1:numel(alpha_ref)
                xline(alpha_ref(ai), 'k:', 'LineWidth', 1.0);
            end
            xlabel('\alpha  (threshold = \alpha \times w)');
            ylabel('Implied N_{iter}  (log scale)');
            title(sprintf('필요 반복 횟수 (경험)  |  order=%d  ε=%.2f  k=%d  [N_{theory}=%d]', ...
                order, e, k_val, N_theory));
            legend('Location', 'northwest');
            hold off;
            exportgraphics(gcf, sprintf('ridge_Nimpl_%s.png', tag), 'Resolution', 200);

            % ---- 수치 요약 출력 ----
            fprintf('\n=== order=%d  eps=%.2f  k=%.1f(=%d)  nT=%d ===\n', ...
                order, e, k_mult, k_val, nT_val);
            fprintf('  w_true=%.3f  w^k(theory)=%.5f  N_theory=%d\n', ...
                w_true, wk_theory, N_theory);
            hdr = sprintf('  %-18s', 'lambda');
            for ai = 1:numel(alpha_ref)
                hdr = [hdr, sprintf('  P(Sj>%.0f%%w)', alpha_ref(ai)*100)]; %#ok<AGROW>
            end
            for ai = 1:numel(alpha_ref)
                hdr = [hdr, sprintf('  P(rc>%.2f)', alpha_ref(ai))]; %#ok<AGROW>
            end
            hdr = [hdr, sprintf('  %10s', 'N(Sj>.8w)')];
            fprintf('%s\n', hdr);
            for li = 1:nLam
                ri_m = idx_all(li);
                row_str = sprintf('  %-18s', lambda_labels{li});
                for ai = 1:numel(alpha_ref)
                    row_str = [row_str, sprintf('  %10.4f', RES(ri_m).P_ref_sj(ai))]; %#ok<AGROW>
                end
                for ai = 1:numel(alpha_ref)
                    row_str = [row_str, sprintf('  %10.4f', RES(ri_m).P_ref_rec(ai))]; %#ok<AGROW>
                end
                p80 = RES(ri_m).P_ref_sj(alpha_ref == 0.80);
                if isempty(p80), p80 = NaN; end
                n80 = ceil(log(1 - conf) / log(max(realmin, 1 - p80)));
                row_str = [row_str, sprintf('  %10d', n80)]; %#ok<AGROW>
                fprintf('%s\n', row_str);
            end
            fprintf('  (이론 N=%d 기준)\n', N_theory);
        end
    end
end
