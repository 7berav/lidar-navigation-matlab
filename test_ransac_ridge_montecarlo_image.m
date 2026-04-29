%% test_ransac_ridge_montecarlo_image.m
% [기존] Ridge vs OLS Monte Carlo 분석 (섹션 3~끝 — 변경 없음)
% [신규] frame_006.ply 단일 피팅 직관 그래프 (PPT용)
%


%% ---- [0] 출력 폴더 / figDir 설정 ----
outDir   = 'test_ridge';
if ~isfolder(outDir), mkdir(outDir); end
runStamp = datestr(now, 'yyyymmdd_HHMM');
figDir   = fullfile(outDir, runStamp);   % test_ridge/20260415_1230/
mkdir(figDir);
fprintf('그림 폴더 → %s\n', figDir);

warning('off', 'MATLAB:nearlySingularMatrix');

%% ---- [1] frame_006 단일 피팅 — lambda=0 (OLS) ----
% 의도:
%   - 차수는 6차 고정 (H6)
%   - 피팅 샘플 개수 = k_mult × nT  (RANSAC이 쓰는 1.4배보다 넉넉한 2배)
%   - lambda = 0 (OLS) 한 번만. lambda 비교는 별도 섹션에서.
% 시각화: figure 각각 따로 (raw / 피팅)
%
% 파라미터 (여기만 수정):
fit_ply          = fullfile('DATA', 'frame_016.ply');
fit_order        = 6;                              % 고정 차수 (H6)
fit_k_mult       = 1.4;% Nfit = k_mult × nT
fit_Nsub_disp    = 3000;                           % scatter3 표시용 subsample 수
fit_view         = [1, 1, 0.2];                   % view 방향 (V3와 동일)
fit_ax_lim       = [-1.2 1.2 -1.2 1.2 -1.2 1.2];   % axis 범위
fit_clim         = [-0.5 0.5];                     % 잔차 컬러 범위
fit_mesh_density = 80;                             % fimplicit3 MeshDensity
fit_surf_color   = [0.90, 0.81, 0.53];             % 곡면 색 (V3와 동일)


% ---- PLY 로드 + 정규화 (load_ply_data.m 방식 동일) ----
%pp_use  = pcread(fit_ply);
%pp_use  = double(pp_use.Location);
%pp_use  = pp_use(all(isfinite(pp_use), 2), :);


totalN = 1500;
PP11 = generateRandomPointsOnHexagonPrism(totalN + 1000) + randn(totalN + 1000, 3) * 0.007;
PP14 = (2 * rand(5000, 3) - 1) * 1.7;

n_out   = round(totalN * 0.1);
n_in    = totalN - n_out;
pp_use = [PP11(1:n_in, :); PP14(1:n_out, :)];




c   = mean(pp_use, 1);
pp_use  = pp_use - c;
%pp_use  = pp_use / 0.1;
% ---- 6차 Fischer 항 구성 + 샘플 개수 결정 ----
TermsF          = homogeneFischerTerms(fit_order);
[FuncsF, ~, ~]  = makeFuncsGradsStack(TermsF);
nTF             = numel(TermsF);
fit_Nsub_fit    = round(fit_k_mult * nTF);         % ≈ 2 × nT

rng(42);
N006        = size(pp_use, 1);
idx_fit     = randperm(N006, min(fit_Nsub_fit,  N006));
PP_fit_sub  = pp_use(idx_fit,  :);    % 피팅용 (≈2×nT 개)
PP_disp     = pp_use;    % scatter3 표시용

% ---- 디자인 행렬 (lambda 공통) ----
A = FuncsF(PP_fit_sub);                  % Nfit × nT
b = ones(size(A, 1), 1);                 % implicit: f = 1

% ---- Figure 1 : raw point cloud ----
figure(1); clf;
scatter3(PP_disp(:,1), PP_disp(:,2), PP_disp(:,3), 3, [0.4 0.6 0.9], 'filled');
axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(fit_view);
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

% ---- Figure 2~ : lambda별 피팅 (각각 독립 figure) ----
fit_lambda_list = [0, 1e-3, 1e-2, 1e-1,3e-1,1,3,10,30,100,300,1000];   % ← 여기서 lambda 조절

for li = 1:numel(fit_lambda_list)
    lam_f  = fit_lambda_list(li);
    [beta_f, ~] = regressionFourthOrder(PP_fit_sub, FuncsF, lam_f, fit_order);   % Sobolev ridge, N-정규화
    f_sym  = TermsF * beta_f;
    Func_f = matlabFunction(f_sym);

    vals_disp = Func_f(PP_disp(:,1), PP_disp(:,2), PP_disp(:,3)) - 1;

    if lam_f == 0
        lam_str = '\lambda=0 (OLS)';
        lam_tag = 'l000';
    else
        lam_str = sprintf('\\lambda=%.0e', lam_f);
        lam_tag = sprintf('l%03d', round(lam_f * 1000));
    end

    figure(1 + li); clf;
    scatter3(PP_disp(:,1), PP_disp(:,2), PP_disp(:,3), 3, vals_disp(:), 'filled');
    hold on;
    fimplicit3(f_sym - 1, ...
        'FaceColor', fit_surf_color, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5, ...
        'MeshDensity', fit_mesh_density);
    hold off;
    colormap(gca, jet); colorbar; clim(fit_clim);
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(lam_str, 'Interpreter', 'tex');
    view(fit_view);
    xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

    exportgraphics(figure(1+li), fullfile(figDir, sprintf('fit_%s_or%d.png', lam_tag, fit_order)), 'Resolution', 200);
end


%% ---- [3] 공통 MC 세팅 ---- (기존 변경 없음)
totalN      = 450;

eps_list    = [0.1];
k_list      = [0.8 1.0 1.2];
%k_list      = [1.0];
order_list  = [6];
Niter       = 2000;
lambda_list  = [1e-3, 1e-2, 1e-1, 3e-1];         % ← 이번 실행에 돌릴 lambda
%lambda_list  = [0, 1e-3,3e-3,1e-2]
%lambda_list  = [ 1e-2,3e-2,1e-1,3e-1]
conf        = 0.95;

% w-스케일 threshold용 alpha 범위
alpha_arr     = linspace(0.10, 1.02, 120);   % thresh = alpha * w_true
alpha_ref_sj  = [0.40, 0.5, 0.60, 0.70];   % preScore bar chart (낮은 range)
alpha_ref_rec = [0.70, 0.80,0.90, 0.95];          % recall bar chart

% 데이터 풀 생성
%PP11 = generateRandomPointsOnHexagonPrism(totalN + 1000) + randn(totalN + 1000, 3) * 0.007;
PP11 = pp_use;
PP14 = (2 * rand(5000, 3) - 1) * 2;

%% ---- [4] RANSAC 파라미터 ---- (기존 변경 없음)
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
    'freezeW',     false,  ...
    'freezeIter',  true,  ...
    'freezeOmega', true   ...
);

%% ---- [5] 실행 루프 ---- (기존 변경 없음)
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
                    lam_tag = regexprep(sprintf('%.0e\', lam), 'e\+0*', 'e');
                    lam_tag = regexprep(lam_tag, 'e-0*', 'm');
                end
                ransacPar.mc.saveVarName = sprintf( ...
                    'MC_ridge_w%02d_k%02d_or%d_%s', ...
                    round(100*e), round(10*k_mult), order, lam_tag);

                fprintf('order=%d  eps=%.2f  k=%.1f(=%d)  lambda=%-8s ... ', ...
                    order, e, k_mult, k_samples, lam_tag);
                t0 = tic;

                [~, beta_best, ~, Log] = PoliNavigationSolver3_3_FischerRansac_MC( ...
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
                RES(row).beta_best  = beta_best;   % 최적 beta (implicit surface용)
            end
        end
    end
end



%% ---- [6] 분석: w-스케일 threshold sweep ---- (기존 변경 없음)
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

%% ---- [7] RES 저장 확인 ---- (기존 변경 없음 — figDir는 [0]에서 생성)
resFile  = fullfile(outDir, sprintf('RES_%s.mat', runStamp));
%save(resFile, 'RES');
fprintf('RES 저장 → %s  (%d 행)\n', resFile, numel(RES));
fprintf('그림 폴더 → %s\n', figDir);

%% ---- [8] 시각화 Fig A / B / C ---- (기존 변경 없음)
% ★ 여기서 lambda_graph 수정 ★
% RES에 없는 값은 자동 skip → 포괄적으로 써두면 됨
lambda_graph = [1e-3, 1e-2, 1e-1, 3e-1, 1e0, 3e0];   % 비교할 lambda 후보 전체

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

            % ---- Fig A: 점수 분포 (preScore + recall 2-panel) ----
            figure(figBase); clf; figBase = figBase + 1;

            subplot(2, 1, 1); hold on; grid on;
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
            xlim([0.3 inf]);
            hold off;

            subplot(2, 1, 2); hold on; grid on;
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
            xlim([0.3 inf]);
            hold off;
            exportgraphics(gcf, fullfile(figDir, sprintf('ridge_dist_%s.png', tag)), 'Resolution', 200);
            savefig(gcf, fullfile(figDir, sprintf('ridge_dist_%s.fig', tag)));

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
                row_str = [row_str, sprintf('  %10d', n_rc80)]; %#ok<AGROW>
                fprintf('%s\n', row_str);
            end
            fprintf('  (이론 N=%d 기준)\n', N_theory);
        end
    end
end

%% ---- [9] MC 최적 beta → scatter3 + fimplicit3 ----
% RANSAC MC에서 나온 best beta로 도형을 직접 그림
% ★ 여기만 바꾸면 됨
target_order  = 4;
target_eps    = 0.00;
target_k_mult = 1.2;
target_lambda_list = [0,1e-3, 3e-3,1e-2];   % 비교할 lambda들 (각각 독립 figure)

% 시각화 파라미터
vis_view        = fit_view;
vis_surf_color  = fit_surf_color;
vis_mesh_density = fit_mesh_density;
vis_clim        = fit_clim;

% PP_vis: 그림에 뿌릴 점군 (MC 데이터 or PLY — 여기서 선택)
% MC 데이터로 그리려면 PPm_use 사용, PLY로 그리려면 pp_use 사용
PP_vis = pp_use;    % ← PLY 점군 기준

for li = 1:numel(target_lambda_list)
    lam_t = target_lambda_list(li);

    hit = find( ...
        [RES.order]  == target_order                   & ...
        abs([RES.eps]    - target_eps)      < 1e-9     & ...
        abs([RES.k_mult] - target_k_mult)   < 1e-9     & ...
        abs([RES.lambda] - lam_t)           < 1e-12);
    if isempty(hit)
        fprintf('[9] 없음: order=%d eps=%.2f k=%.1f lambda=%.0e → skip\n', ...
            target_order, target_eps, target_k_mult, lam_t);
        continue;
    end
    ri = hit(end);

    beta9  = RES(ri).beta_best;
    Terms9 = homogeneFischerTerms(target_order);
    f_sym9 = Terms9 * beta9(:);
    Func9  = matlabFunction(f_sym9);

    vals9 = Func9(PP_vis(:,1), PP_vis(:,2), PP_vis(:,3)) - 1;

    if lam_t == 0
        lam_str = '\lambda=0 (OLS)';
    else
        lam_str = sprintf('\\lambda=%.0e', lam_t);
    end

    figure; clf;
    scatter3(PP_vis(:,1), PP_vis(:,2), PP_vis(:,3), 3, vals9(:), 'filled');
    hold on;
    fimplicit3(f_sym9 - 1, ...
        'FaceColor', vis_surf_color, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5, ...
        'MeshDensity', vis_mesh_density);
    hold off;
    colormap(gca, jet); colorbar; caxis(vis_clim);
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('MC best beta  |  H%d  %s', target_order, lam_str), 'Interpreter', 'tex');
    view(vis_view);
end
