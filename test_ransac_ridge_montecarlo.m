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

eps_list    = [0.1, 0.20];
k_list      = [1.0 , 1.2, 1.4, 1.6, 1.8];
%k_list      = [1.0];
order_list  = [ 4];
Niter       = 10000;
lambda_list  = [0, 1e-3, 2e-3, 3e-3, 1e-2, 3e-2, 1e-1];         % ← 이번 실행에 돌릴 lambda

conf        = 0.95;

% w-스케일 threshold용 alpha 범위
alpha_arr     = linspace(0.10, 1.02, 120);   % thresh = alpha * w_true
alpha_ref_sj  = [0.40, 0.5, 0.60, 0.70];   % preScore bar chart (낮은 range)
alpha_ref_rec = [0.70, 0.80,0.90, 0.95];          % recall bar chart

% 데이터 풀 생성
PP11 = generateRandomPointsOnHexagonPrism(totalN + 1000) + randn(totalN + 1000, 3) * 0.007;
PP14 = (2 * rand(5000, 3) - 1) * 2;

%% ---- RANSAC 파라미터 ----
ransacPar = struct( ...
    'maxIter',        Niter,  ...
    'conf',           conf,   ...
    'thresh',         0.45,   ...
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
% RES는 누적 방식 → 섹션 재실행해도 기존 결과 유지
% 완전 초기화가 필요하면 워크스페이스에서 직접 clear RES
if ~exist('RES', 'var') || isempty(RES)
    RES = struct([]);
end
row = numel(RES);   % 기존 항목 이어서 추가
%%
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
                    lam_tag = regexprep(sprintf('%.0e', lam), 'e\+0*', 'e');
                    lam_tag = regexprep(lam_tag, 'e-0*', 'm');
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
                mseAlls    = double([Log.mseAll]');
                mseInliers = double([Log.mseInlier]');

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
                RES(row).mseAlls    = mseAlls;
                RES(row).mseInliers = mseInliers;
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

    % alpha_ref 기준 상위 비율 (bar chart용, sj/rec 별도 range)
    P_ref_sj  = arrayfun(@(a) mean(valid_ps > a * w),  alpha_ref_sj);
    P_ref_rec = arrayfun(@(a) mean(valid_rc > a),       alpha_ref_rec);

    RES(ri).alpha_arr  = alpha_arr;
    RES(ri).thresh_sj  = thresh_sj;
    RES(ri).thresh_rec = thresh_rec;
    RES(ri).P_sj       = P_sj;
    RES(ri).P_rc       = P_rc;
    RES(ri).N_impl_sj  = N_impl_sj;
    RES(ri).P_ref_sj   = P_ref_sj;   % (1 x numel(alpha_ref_sj))
    RES(ri).P_ref_rec  = P_ref_rec;  % (1 x numel(alpha_ref_rec))
end

%% ---- RES 저장 + 출력 폴더 ----
outDir = 'test_ridge';
if ~isfolder(outDir), mkdir(outDir); end
runStamp = datestr(now, 'yyyymmdd_HHMM');
figDir   = fullfile(outDir, runStamp);       % test_ridge/20260324_1530/
mkdir(figDir);
resFile  = fullfile(outDir, sprintf('RES_%s.mat', runStamp));
%save(resFile, 'RES');
fprintf('RES 저장 → %s  (%d 행)\n', resFile, numel(RES));
fprintf('그림 폴더 → %s\n', figDir);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            %% ---- 시각화 ----
% ★ 여기서 lambda_graph 수정 ★
% RES에 없는 값은 자동 skip → 포괄적으로 써두면 됨
lambda_graph = [0, 1e-3, 2e-3, 3e-3, 1e-2, 3e-2];   % 비교할 lambda 후보 전체

% RES에서 실제 존재하는 조합 추출 → 공통세팅 의존 제거
orders_graph = unique([RES.order]);
eps_graph    = unique([RES.eps]);
k_graph      = unique([RES.k_mult]);

colors = lines(numel(lambda_graph));
lambda_labels = cell(1, numel(lambda_graph));
for li = 1:numel(lambda_graph)
    if lambda_graph(li) == 0
        lambda_labels{li} = 'OLS (\lambda=0)';
    else
        lambda_labels{li} = sprintf('\\lambda=%.0e', lambda_graph(li));
    end
end

figBase = 20;
% Fig A y축 통일용 컨테이너 (모든 세팅 그린 후 일괄 적용)
ax_figA_ps  = {};
ax_figA_rc  = {};
ax_figA_mse = {};
fig_A_handles = {};
fig_A_tags    = {};

for order = orders_graph
    for e = eps_graph
        for k_mult = k_graph

            % ── 이 (order, eps, k_mult) 조합에 대해 lambda별로 RES 행 직접 찾기 ──
            % lambda까지 명시적으로 매칭 → 순서 가정 없음, RES가 누적돼도 안전
            % lambda_graph 에 있는 lambda별로 RES 행 직접 찾기
            % 없는 lambda는 skip → lambda_graph 일부만 RES에 있어도 동작
            nLam       = numel(lambda_graph);
            ri_per_lam = zeros(1, nLam);
            keep       = true(1, nLam);
            for li = 1:nLam
                hit = find( ...
                    [RES.order]  == order                         & ...
                    abs([RES.eps]    - e)                < 1e-9   & ...
                    abs([RES.k_mult] - k_mult)           < 1e-9   & ...
                    abs([RES.lambda] - lambda_graph(li)) < 1e-12);
                if isempty(hit)
                    fprintf('  [없음] lambda=%.0e → skip\n', lambda_graph(li));
                    keep(li) = false;
                else
                    ri_per_lam(li) = hit(end);
                end
            end
            % 유효한 lambda만 추출
            ri_per_lam  = ri_per_lam(keep);
            lam_colors  = colors(keep, :);
            lam_labels  = lambda_labels(keep);
            nLam        = sum(keep);
            if nLam == 0, continue; end

            % 공통 메타데이터
            ri_ref    = ri_per_lam(1);
            nT_val    = RES(ri_ref).nT;
            k_val     = RES(ri_ref).k;
            w_true    = RES(ri_ref).w_true;
            wk_theory = RES(ri_ref).wk_theory;
            N_theory  = RES(ri_ref).N_theory;
            tag       = sprintf('or%d_e%02d_k%02d', order, round(100*e), round(10*k_mult));

            % ---- Fig A: 점수 분포 (preScore + recall + mseInlier 3-panel) ----
            % 저장은 y축 통일 후 일괄 처리 (루프 끝 참고)
            fig_A = figure(figBase); clf; figBase = figBase + 1;

            ax_ps = subplot(3, 1, 1); hold on; grid on;
            for li = 1:nLam
                ps = RES(ri_per_lam(li)).preScores;
                histogram(ps(isfinite(ps)), 60, ...
                    'FaceAlpha', 0.35, 'EdgeAlpha', 0.15, ...
                    'FaceColor', lam_colors(li, :), 'DisplayName', lam_labels{li});
            end
            for ai = 1:numel(alpha_ref_sj)
                xline(alpha_ref_sj(ai) * w_true, 'k:', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%.0f%%×w', alpha_ref_sj(ai)*100));
            end
            xline(w_true, 'r-', 'LineWidth', 1.5, 'DisplayName', 'w (이론 최대)');
            xlabel('preScore (Sj)'); ylabel('count');
            title(sprintf('preScore 분포  |  order=%d  ε=%.2f  k=%.1f×%d=%d', ...
                order, e, k_mult, nT_val, k_val));
            legend('Location', 'northeast', 'FontSize', 7);
            xlim([0.1 0.4]);
            hold off;
            ax_figA_ps{end+1} = ax_ps;

            ax_rc = subplot(3, 1, 2); hold on; grid on;
            for li = 1:nLam
                rc = RES(ri_per_lam(li)).recalls;
                histogram(rc(isfinite(rc)), 60, ...
                    'FaceAlpha', 0.35, 'EdgeAlpha', 0.15, ...
                    'FaceColor', lam_colors(li, :), 'DisplayName', lam_labels{li});
            end
            for ai = 1:numel(alpha_ref_rec)
                xline(alpha_ref_rec(ai), 'k:', 'LineWidth', 1.2);
            end
            xlabel('recall'); ylabel('count');
            title('recall 분포  |  updateThresh=0 → 전 iter 로컬 opt');
            legend('Location', 'northwest', 'FontSize', 7);
            xlim([0.3 1]);
            hold off;
            ax_figA_rc{end+1} = ax_rc;

            ax_mse = subplot(3, 1, 3); hold on; grid on;
            for li = 1:nLam
                ms = RES(ri_per_lam(li)).mseInliers;
                histogram(ms(isfinite(ms)), 60, ...
                    'FaceAlpha', 0.35, 'EdgeAlpha', 0.15, ...
                    'FaceColor', lam_colors(li, :), 'DisplayName', lam_labels{li});
            end
            xlabel('mseInlier  (mean r^2, 인라이어)'); ylabel('count');
            title('inlier MSE 분포');
            legend('Location', 'northeast', 'FontSize', 7);
            xlim([0.02 0.08]);
            hold off;
            ax_figA_mse{end+1} = ax_mse;

            fig_A_handles{end+1} = fig_A;
            fig_A_tags{end+1}    = tag;

            % ---- Fig B: 상위 비율 bar chart ----
            figure(figBase); clf; figBase = figBase + 1;

            barData_sj  = zeros(numel(alpha_ref_sj),  nLam);
            barData_rec = zeros(numel(alpha_ref_rec), nLam);
            for li = 1:nLam
                barData_sj(:,  li) = RES(ri_per_lam(li)).P_ref_sj';
                barData_rec(:, li) = RES(ri_per_lam(li)).P_ref_rec';
            end

            subplot(1, 2, 1);
            b = bar(barData_sj);
            for li = 1:nLam
                b(li).FaceColor = lam_colors(li, :);
                b(li).DisplayName = lam_labels{li};
            end
            set(gca, 'XTickLabel', arrayfun(@(a) sprintf('%.0f%%×w', a*100), ...
                alpha_ref_sj, 'UniformOutput', false));
            ylabel('P(Sj > threshold)'); title('preScore 상위 비율');
            legend('Location', 'northeast', 'FontSize', 7);
            grid on; ylim([0 min(1, max(barData_sj(:))*1.25 + 0.02)]);

            subplot(1, 2, 2);
            b = bar(barData_rec);
            for li = 1:nLam
                b(li).FaceColor = lam_colors(li, :);
                b(li).DisplayName = lam_labels{li};
            end
            set(gca, 'XTickLabel', arrayfun(@(a) sprintf('recall>%.2f', a), ...
                alpha_ref_rec, 'UniformOutput', false));
            ylabel('P(recall > threshold)'); title('recall 상위 비율');
            legend('Location', 'northeast', 'FontSize', 7);
            grid on; ylim([0 min(1, max(barData_rec(:))*1.25 + 0.02)]);
            sgtitle(sprintf('상위 비율  |  ε=%.2f  k=%d  nT=%d', e, k_val, nT_val));
            exportgraphics(gcf, fullfile(figDir, sprintf('ridge_toprate_%s.png', tag)), 'Resolution', 200);
            savefig(gcf, fullfile(figDir, sprintf('ridge_toprate_%s.fig', tag)));

            % ---- Fig C: implied N_iter vs alpha ----
            figure(figBase); clf; figBase = figBase + 1;
            hold on; grid on;
            for li = 1:nLam
                semilogy(RES(ri_per_lam(li)).alpha_arr, RES(ri_per_lam(li)).N_impl_sj, '-', ...
                    'Color', lam_colors(li, :), 'LineWidth', 1.8, 'DisplayName', lam_labels{li});
            end
            for ai = 1:numel(alpha_ref_sj)
                xline(alpha_ref_sj(ai), 'k:', 'LineWidth', 1.0);
            end
            xlabel('\alpha  (threshold = \alpha \times w)');
            ylabel('Implied N_{iter}  (log scale)');
            title(sprintf('필요 반복 횟수 (경험)  |  order=%d  ε=%.2f  k=%d  [N_{theory}=%d]', ...
                order, e, k_val, N_theory));
            legend('Location', 'northwest');
            hold off;
            exportgraphics(gcf, fullfile(figDir, sprintf('ridge_Nimpl_%s.png', tag)), 'Resolution', 200);
            savefig(gcf, fullfile(figDir, sprintf('ridge_Nimpl_%s.fig', tag)));

            % ---- 수치 요약 출력 ----
            fprintf('\n=== order=%d  eps=%.2f  k=%.1f(=%d)  nT=%d ===\n', ...
                order, e, k_mult, k_val, nT_val);
            fprintf('  w_true=%.3f  w^k(theory)=%.5f  N_theory=%d\n', ...
                w_true, wk_theory, N_theory);
            hdr = sprintf('  %-18s', 'lambda');
            for ai = 1:numel(alpha_ref_sj)
                hdr = [hdr, sprintf('  P(Sj>%.0f%%w)', alpha_ref_sj(ai)*100)]; %#ok<AGROW>
            end
            for ai = 1:numel(alpha_ref_rec)
                hdr = [hdr, sprintf('  P(rc>%.2f)', alpha_ref_rec(ai))]; %#ok<AGROW>
            end
            hdr = [hdr, sprintf('  %10s', 'N(rc>.80)')];
            fprintf('%s\n', hdr);
            for li = 1:nLam
                ri_m = ri_per_lam(li);
                row_str = sprintf('  %-18s', lam_labels{li});
                for ai = 1:numel(alpha_ref_sj)
                    row_str = [row_str, sprintf('  %10.4f', RES(ri_m).P_ref_sj(ai))]; %#ok<AGROW>
                end
                for ai = 1:numel(alpha_ref_rec)
                    row_str = [row_str, sprintf('  %10.4f', RES(ri_m).P_ref_rec(ai))]; %#ok<AGROW>
                end
                p_rc80 = RES(ri_m).P_ref_rec(alpha_ref_rec == 0.80);
                if isempty(p_rc80), p_rc80 = NaN; end
                n_rc80 = ceil(log(1 - conf) / log(max(realmin, 1 - p_rc80)));
                row_str = [row_str, sprintf('  %10d', n_rc80)]; %#ok!giAGROW>
                fprintf('%s\n', row_str);
            end
            fprintf('  (이론 N=%d 기준)\n', N_theory);
        end
    end
end

%% ---- Fig A y축 통일 + 저장 ----
% 세 패널(preScore / recall / mseInlier) 모두 동일한 ylim 적용
% MATLAB 자동 패딩 제거: histogram BinCounts 직접 추출 → 실제 피크 * 1.05
if ~isempty(ax_figA_ps)
    % 실제 bin 최대값 추출 (auto-padding 없이)
    histMax = @(ax) max(cellfun(@(h) max([h.BinCounts(:); 0]), ...
        num2cell(findobj(ax, 'Type', 'histogram'))));

    raw_ps  = max(cellfun(histMax, ax_figA_ps));
    raw_rc  = max(cellfun(histMax, ax_figA_rc));
    raw_mse = max(cellfun(histMax, ax_figA_mse));

    % 세 패널 전체를 단일 ylim으로 통일 (5% 여유)
    ymax_all = ceil(max([raw_ps, raw_rc, raw_mse]) * 1.05);

    for i = 1:numel(fig_A_handles)
        ylim(ax_figA_ps{i},  [0 ymax_all]);
        ylim(ax_figA_rc{i},  [0 ymax_all]);
        ylim(ax_figA_mse{i}, [0 ymax_all]);
        exportgraphics(fig_A_handles{i}, ...
            fullfile(figDir, sprintf('ridge_dist_%s.png', fig_A_tags{i})), 'Resolution', 200);
        savefig(fig_A_handles{i}, ...
            fullfile(figDir, sprintf('ridge_dist_%s.fig', fig_A_tags{i})));
    end
end
