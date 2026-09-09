%% aggregate_shapes.m — 도형별 lcurve_demo 결과 통합 분석 (Sec 3.1 마무리)
%
% 각 도형의 최신 test_ridge/<stamp>_<shape>/lcurve_demo_results.csv 를 모아:
%   표 1   : 도형×occ별 λ*(CD_all-opt), df*_sub(뷰 중앙값), df*_full, 괴리
%   표 2   : GCV-선택 λ의 CD regret (subset, GT-free 규칙 유효성)
%   그림 A : df* 괴리 vs occlusion, 도형 오버레이 (크기 불변성 붕괴점)
%   그림 B : CD_all vs df(λ 파라메트릭), 도형 오버레이 — df가 기하 품질 지배
%   저장   : test_ridge/sec31_interface.mat (3.2 인터페이스: df*, 괴리, σ̂)
%
% 게이트 판정(콘솔): (a) λ* 내부점 여부  (b) GCV regret ≤ 0.02  (c) df* basin 겹침

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))

shapes  = {'hexagon', 'cube', 'cylinder', 'sphere', 'ellipsoid'};
kFocus  = 1.4;
outRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'test_ridge');

aggDir = fullfile(outRoot, sprintf('aggregate_%s', datestr(now, 'yyyymmdd_HHMM')));
mkdir(aggDir);

nS = numel(shapes);
S  = struct([]);   % 3.2 인터페이스 컨테이너

for is = 1:nS
    shape = shapes{is};
    % 최신 폴더 탐색 (타임스탬프 사전순 최대 = 최신, _quick 제외)
    dd   = dir(outRoot);
    nm   = {dd([dd.isdir]).name};
    pat  = sprintf('^\\d{8}_\\d{4}_%s(_or\\d+)?$', shape);
    hits = nm(~cellfun('isempty', regexp(nm, pat, 'once')));
    assert(~isempty(hits), '결과 폴더 없음: %s (run_all_shapes 먼저 실행)', shape);
    hits   = sort(hits);
    folder = fullfile(outRoot, hits{end});
    T = readtable(fullfile(folder, 'lcurve_demo_results.csv'));
    fprintf('[%s] %s  (%d rows)\n', shape, hits{end}, height(T));

    occs = unique(T.occ)';
    lams = unique(T.lambda)';
    nO = numel(occs); nL = numel(lams);
    views = unique(T.view)'; nV = numel(views);

    lamStarSub = nan(1,nO); dfStarMed = nan(1,nO); dfStarIQR = nan(1,nO);
    lamStarFul = nan(1,nO); dfFulStar = nan(1,nO);
    regretMed  = nan(1,nO); regretMax = nan(1,nO);
    sigmaHat   = nan(1,nO);
    % 그림 B용 λ-파라메트릭 곡선 (뷰 평균)
    dfCurveSub = nan(nO,nL); cdCurveSub = nan(nO,nL);
    dfCurveFul = nan(nO,nL); cdCurveFul = nan(nO,nL);

    for io = 1:nO
        occ = occs(io);
        mSub = strcmp(T.fit_type,'subset') & abs(T.occ-occ)<1e-9 & abs(T.k-kFocus)<1e-9;
        mFul = strcmp(T.fit_type,'full')   & abs(T.occ-occ)<1e-9;

        % λ-파라메트릭 뷰 평균 곡선
        for il = 1:nL
            ms = mSub & abs(T.lambda-lams(il))<1e-15;
            mf = mFul & abs(T.lambda-lams(il))<1e-15;
            dfCurveSub(io,il) = mean(T.df(ms),'omitnan');
            cdCurveSub(io,il) = mean(T.cd_all(ms),'omitnan');
            dfCurveFul(io,il) = mean(T.df(mf),'omitnan');
            cdCurveFul(io,il) = mean(T.cd_all(mf),'omitnan');
        end

        % λ* : 뷰 평균 CD_all의 argmin
        [~,iaS] = min(cdCurveSub(io,:)); lamStarSub(io) = lams(iaS);
        [~,iaF] = min(cdCurveFul(io,:)); lamStarFul(io) = lams(iaF);
        dfFulStar(io) = dfCurveFul(io,iaF);

        % df* 앙상블: 뷰별 CD-opt에서의 df → median/IQR (docs/20260721 §3.1)
        dfV = nan(1,nV); regV = nan(1,nV); mseV = nan(1,nV);
        for iv = 1:nV
            tt = sortrows(T(mSub & T.view==views(iv), :), 'lambda');
            if isempty(tt), continue; end
            [~,ic] = min(tt.cd_all);
            dfV(iv) = tt.df(ic);
            [~,ig] = min(tt.gcv);
            regV(iv) = tt.cd_all(ig) - tt.cd_all(ic);   % GCV의 CD regret
            mseV(iv) = tt.mse(ic);                       % λ*에서의 대수 MSE
        end
        dfStarMed(io) = median(dfV,'omitnan');
        dfStarIQR(io) = iqr(dfV);
        regretMed(io) = median(regV,'omitnan');
        regretMax(io) = max(regV,[],'omitnan');
        % σ̂: λ*에서의 subset 대수 잔차 스케일 프록시 (뷰 중앙값 MSE의 제곱근).
        % raw residual이 CSV에 없어 MAD 기반은 3.2 solver 내부에서 재계산 예정.
        sigmaHat(io) = sqrt(median(mseV,'omitnan'));
    end

    S(is).shape       = shape;
    S(is).folder      = folder;
    S(is).occ_list    = occs;
    S(is).lambda_list = lams;
    S(is).lambda_star_sub = lamStarSub;
    S(is).lambda_star_full= lamStarFul;
    S(is).df_star_med = dfStarMed;
    S(is).df_star_iqr = dfStarIQR;
    S(is).df_full_star= dfFulStar;
    S(is).gap         = dfStarMed - dfFulStar;
    S(is).sigma_hat   = sigmaHat;
    S(is).regret_med  = regretMed;
    S(is).regret_max  = regretMax;
    S(is).df_curve_sub = dfCurveSub;  S(is).cd_curve_sub = cdCurveSub;
    S(is).df_curve_ful = dfCurveFul;  S(is).cd_curve_ful = cdCurveFul;
end

occs = S(1).occ_list;  lams = S(1).lambda_list;  nO = numel(occs);

%% ---- 표 1: λ*, df*, 괴리 ----
fprintf('\n===== 표 1. 도형×occ: λ*(CD_all-opt), df*_sub(med±IQR), df*_full, 괴리 =====\n');
fprintf('%-10s %-5s | %-9s %-14s %-9s %-7s\n', 'shape','occ','λ*_sub','df*_sub','df*_full','gap');
for is = 1:nS
    for io = 1:nO
        fprintf('%-10s %-5.2f | %-9.0e %6.2f ± %-5.2f %-9.2f %+-7.2f\n', ...
            S(is).shape, occs(io), S(is).lambda_star_sub(io), ...
            S(is).df_star_med(io), S(is).df_star_iqr(io), ...
            S(is).df_full_star(io), S(is).gap(io));
    end
end

%% ---- 표 2: GCV regret ----
fprintf('\n===== 표 2. GCV-선택 λ의 CD_all regret (subset k=%.1f) =====\n', kFocus);
fprintf('%-10s', 'shape');
fprintf(' occ=%-11.2f', occs); fprintf('\n');
for is = 1:nS
    fprintf('%-10s', S(is).shape);
    for io = 1:nO
        fprintf(' %+.3f/%+.3f  ', S(is).regret_med(io), S(is).regret_max(io));
    end
    fprintf('\n');
end
fprintf('(각 셀: median/max over views)\n');

%% ---- 그림 A: df* 괴리 vs occlusion ----
cols = lines(nS);
fig = figure('Visible','off','Position',[60 60 760 540]); hold on; grid on;
for is = 1:nS
    plot(100*occs, S(is).gap, 'o-', 'Color', cols(is,:), ...
        'LineWidth', 1.8, 'MarkerFaceColor', cols(is,:), 'DisplayName', S(is).shape);
end
yline(0, 'k:', 'HandleVisibility','off');
xlabel('occlusion ratio [%]');
ylabel('df^*_{sub} - df^*_{full}  (CD_{all}-opt 기준)');
title(sprintf('크기 불변성 진단: df^* 괴리 vs occlusion  (subset k=%.1f)', kFocus));
legend('Location','best');
exportgraphics(fig, fullfile(aggDir,'aggA_df_gap_vs_occ.png'), 'Resolution', 200);
savefig(fig, fullfile(aggDir,'aggA_df_gap_vs_occ.fig'));

%% ---- 그림 B: CD_all vs df (λ 파라메트릭) ----
occShow = intersect([0.30 0.50 0.70], occs, 'stable');
fig = figure('Visible','off','Position',[40 40 1400 460]);
for p = 1:numel(occShow)
    io = find(abs(occs - occShow(p)) < 1e-9, 1);
    subplot(1, numel(occShow), p); hold on; grid on;
    for is = 1:nS
        plot(S(is).df_curve_sub(io,:), S(is).cd_curve_sub(io,:), 'o-', ...
            'Color', cols(is,:), 'LineWidth', 1.5, 'DisplayName', S(is).shape);
        plot(S(is).df_curve_ful(io,:), S(is).cd_curve_ful(io,:), '--', ...
            'Color', cols(is,:), 'LineWidth', 1.2, 'HandleVisibility','off');
        [~,im] = min(S(is).cd_curve_sub(io,:));
        plot(S(is).df_curve_sub(io,im), S(is).cd_curve_sub(io,im), 'p', ...
            'Color', cols(is,:), 'MarkerFaceColor', cols(is,:), ...
            'MarkerSize', 11, 'HandleVisibility','off');
    end
    xlabel('df(\lambda) = tr(H_\lambda)'); ylabel('CD_{all}');
    title(sprintf('occ=%d%%', round(100*occShow(p))));
    if p == 1, legend('Location','best','FontSize',7); end
end
sgtitle(sprintf(['CD_{all} vs effective dof — solid: subset k=%.1f, ' ...
    'dashed: all visible, star: CD-opt'], kFocus));
exportgraphics(fig, fullfile(aggDir,'aggB_cd_vs_df.png'), 'Resolution', 200);
savefig(fig, fullfile(aggDir,'aggB_cd_vs_df.fig'));

%% ---- 3.2 인터페이스 저장 ----
ifaceFile = fullfile(outRoot, 'sec31_interface.mat');
save(ifaceFile, 'S', 'kFocus');
fprintf('\n3.2 인터페이스 저장 → %s\n', ifaceFile);

%% ---- 게이트 판정 ----
fprintf('\n===== 3.2 진입 게이트 =====\n');
lamMax = max(lams);
edgeHit = false; regretBad = false;
for is = 1:nS
    onEdge = abs(S(is).lambda_star_sub - lamMax) < 1e-15;
    if any(onEdge)
        edgeHit = true;
        fprintf('  [경계] %s: occ=%s 에서 λ*가 그리드 최대(%.0e)\n', S(is).shape, ...
            mat2str(round(100*occs(onEdge))), lamMax);
    end
    bad = S(is).regret_med > 0.02;
    if any(bad)
        regretBad = true;
        fprintf('  [regret] %s: occ=%s 에서 GCV regret(med) > 0.02\n', S(is).shape, ...
            mat2str(round(100*occs(bad))));
    end
end
if ~edgeHit,   fprintf('  (a) λ* 전부 내부점 — 그리드 충분 ✓\n'); end
if ~regretBad, fprintf('  (b) GCV regret(med) ≤ 0.02 전 도형·occ ✓\n'); end
% (c) df* basin 겹침: occ별 도형 간 df* 범위
fprintf('  (c) df* basin (occ별 도형 간 [min..max]):\n');
for io = 1:nO
    v = arrayfun(@(s) s.df_star_med(io), S);
    fprintf('      occ=%2d%%: [%5.2f .. %5.2f]  spread=%.2f\n', ...
        round(100*occs(io)), min(v), max(v), max(v)-min(v));
end
fprintf('\n통합 그림/표 저장 → %s\n', aggDir);
