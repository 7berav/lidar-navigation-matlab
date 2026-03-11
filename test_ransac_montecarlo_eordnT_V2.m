%% 공통 세팅
totalN       = 10000;
eps_list     = [0.0 0.10 0.20 0.30]; %[0.00 0.10 0.20 0.30];
k_list       = [1.2 1.4];            % nT 배수 (저장명에는 round(10*k) 사용)
order_list   = [ 4 ];
Niter        = 2000;

PP11 = generateRandomPointsOnHexagonPrism(10400) + randn(10400,3)*0.007;  % inlier pool
PP14 = (2*rand(5000,3)-1)*2;                                             % outlier pool (Uniform[-2,2])
PP12 = generateRandomPointsOnCylinder(10400) + randn(10400,3)*0.007;  % inlier pool
PP13 = generateRandomPointsOnCube(10400) + randn(10400,3)*0.007;  % inlier pool

% 고정 RANSAC 파라미터(필요 필드만 세팅)
ransacPar = struct( ...
  'maxIter', Niter, ...
  'conf', 0.96, ...
  'thresh', 0.40, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0, ...
  'locIters', 4, ...
  'momentum', 0.66, ...
  'damping', 0.85, ...
  'reg', 1e-6, ...
  'k', 1.4 ...              % 루프에서 덮어씀
);
ransacPar.mc  = struct('on', true, 'saveVarName', 'MC_scores_tmp', 'saveMatFile','', 'time', true, 'localOff', false);
ransacPar.sim = struct('freezeW', false, 'freezeIter', true, 'freezeOmega', false);

warning('off','MATLAB:nearlySingularMatrix');
%%
ResRAN = struct([]);
row    = 0;

%%

for order = order_list
    % ----- order별 basis/미분 생성 -----
    TermsC = homogeneFischerTerms(order);
    [Funcs1, Grads1, ~] = makeFuncsGradsStack(TermsC);
    nT = numel(TermsC);   
    for e = eps_list
        n_out = round(totalN*e);  n_in = totalN - n_out;
        if n_in > size(PP13,1) || n_out > size(PP14,1)
            warning('데이터 풀 부족: order=%d, e=%.2f (n_in=%d, n_out=%d) → 건너뜀', order, e, n_in, n_out);
            continue;
        end
        PPm_use = [ PP13(1:n_in,:); PP14(1:n_out,:) ];

        for k = k_list
            % ---- 필드만 수정 ----
            %ransacPar.maxIter = Niter;
            ransacPar.k       = k;

            % 저장 변수명: MC_scores_w%02d_nT%d_N%d_or%d
            ransacPar.mc.saveVarName = sprintf('MC_scores_w%02d_nT%d_N%d_or%d_cub', ...
                                               round(100*e), round(10*k), Niter, order);

            % 호출 (DisplacementLocal 내부에서 ransacPar 전달되도록 되어 있어야 함)
            [DispRAN, BetaRAN, inlierMaskRAN,Log] = ...
                PoliNavigationSolver3_3_FischerRansac_MC( ...
                    0, PPm_use, order, numel(TermsC), Funcs1, Grads1, ransacPar); 
                            % 외부 SDF (signed)
            P_result = PPm_use - DispRAN;   % 모델 좌표계
            Phi      = Funcs1([P_result(:,1), P_result(:,2), P_result(:,3)]);  % N×nT
            Fval     = Phi * BetaRAN - 1;                                       % N×1
            Gstack   = Grads1(P_result);                                     % (3N)×nT
            gv       = Gstack * BetaRAN;  Nn = size(P_result,1);
            gn       = sqrt( max(1e-12, gv(1:Nn).^2 + gv(Nn+1:2*Nn).^2 + gv(2*Nn+1:3*Nn).^2) );
            d_sdf    = Fval ./ gn;
            
            inlierMask_all = abs(Fval) <= ransacPar.thresh;
            SDF_RMS  = sqrt(mean(d_sdf.^2));
            SDF_P50  = prctile(d_sdf, 50);



        end
    end
end
%%
%{

for order = order_list
    % ----- order별 basis/미분 생성 -----
    TermsC = homogeneFischerTerms(order);
    [Funcs1, Grads1, ~] = makeFuncsGradsStack(TermsC);
    nT = numel(TermsC);   
    for e = eps_list
        n_out = round(totalN*e);  n_in = totalN - n_out;
        if n_in > size(PP12,1) || n_out > size(PP14,1)
            warning('데이터 풀 부족: order=%d, e=%.2f (n_in=%d, n_out=%d) → 건너뜀', order, e, n_in, n_out);
            continue;
        end
        PPm_use = [ PP12(1:n_in,:); PP14(1:n_out,:) ];

        for k = k_list
            % ---- 필드만 수정 ----
            ransacPar.maxIter = Niter;
            ransacPar.k       = k;

            % 저장 변수명: MC_scores_w%02d_nT%d_N%d_or%d
            ransacPar.mc.saveVarName = sprintf('MC_scores_w%02d_nT%d_N%d_or%d_cyl', ...
                                               round(100*e), round(10*k), Niter, order);

            % 호출 (DisplacementLocal 내부에서 ransacPar 전달되도록 되어 있어야 함)
            [DispRAN, BetaRAN, inlierMaskRAN] = ...
                PoliNavigationSolver3_2_FischerRansac_MC( ...
                    0, PPm_use, order, numel(TermsC), Funcs1, Grads1, ransacPar); 
        end
    end
end
%}
%%

order_list  = [ 4 ];                    % ← order 먼저
w_list_pct  = [0 10 20 30];
k_list_pct  = [12 14];
Niter       = 2000;
%timelim     = [10 20 60];
timelim     = [20];
outdir = fullfile(pwd, 'image_cub/enT_246_4');
if ~exist(outdir,'dir'), mkdir(outdir); end

for oi = 1:numel(order_list)
  ordr = order_list(oi);

  for wi = 1:numel(w_list_pct)
    for ki = 1:numel(k_list_pct)
      wpct = w_list_pct(wi);
      kpct = k_list_pct(ki);

      % 변수명/태그: MC_scores_w%02d_nT%d_N%d_or%d
      mcvar = sprintf('MC_scores_w%02d_nT%d_N%d_or%d_cub', wpct, kpct, Niter, ordr);
      tag   = erase(mcvar, 'MC_scores_');

      if ~evalin('base', sprintf('exist(''%s'',''var'')', mcvar))
        warning('없음: %s (건너뜀)', mcvar);
        continue;
      end
      X = evalin('base', mcvar);   % [Sj, SjLoc, raw_ms, loc_ms] (3,4열 없을 수 있음)

      % ===== Fig1: score 히스토그램 =====
      scores_raw = X(:,1); m1 = isfinite(scores_raw);
      scores_loc = X(:,2); m2 = isfinite(scores_loc);

      if any(m1|m2)
        [~, edges] = histcounts([scores_raw(m1); scores_loc(m2)], 'BinMethod','fd');
        p95_score_raw = prctile(scores_raw(m1), 95);
        p95_score_loc = NaN; if any(m2), p95_score_loc = prctile(scores_loc(m2), 95); end

        figure(1); clf;
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        nexttile;
        histogram(scores_raw(m1), edges); grid on; title(sprintf('Sj (raw) | or=%d', ordr));
        xlabel('score'); ylabel('count');
        hold on;
        xline(p95_score_raw, 'r', 'LineWidth', 1.5);
        yl = ylim; text(min(max(p95_score_raw,0),0.7), yl(2)*0.9, sprintf('95%% = %.4g', p95_score_raw), ...
                        'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
        xlim([0 0.7]); hold off;

        nexttile;
        histogram(scores_loc(m2), edges); grid on; title('SjLoc (local)');
        xlabel('score'); ylabel('count');
        hold on;
        if isfinite(p95_score_loc)
          xline(p95_score_loc, 'r', 'LineWidth', 1.5);
          yl = ylim;
          xpos = min(max(p95_score_loc, 0), 0.7);
          text(xpos, yl(2)*0.9, sprintf('95%% = %.4g', p95_score_loc), ...
               'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
        end
        xlim([0 0.7]); hold off;

        exportgraphics(figure(1), fullfile(outdir, sprintf('score_%s.png', tag)), 'Resolution', 300);
        savefig(1,              fullfile(outdir, sprintf('score_%s.fig', tag)));
      end

      % ===== Fig2: Sj vs SjLoc 산점도 (원래 그리기 그대로, 선형축) =====
      mask = m1 & m2;
      topP = 90; topraw = 99.8; thresh = 0.6;
      if any(mask)
        xr = prctile(scores_raw(mask), [topraw 100]); xr = thresh * xr;
        yr = prctile(scores_loc(mask), [topP   100]);

        zoomMask   = mask & scores_loc >= yr(1);
        threshMask = mask & scores_raw >= xr(1);

        figure(2); clf; hold on; grid on;
        plot(scores_raw(mask),   scores_loc(mask),   '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6);
        plot(scores_raw(zoomMask),   scores_loc(zoomMask),   '.', 'MarkerSize', 10);
        plot(scores_raw(threshMask), scores_loc(threshMask), '.', 'MarkerSize', 5);
        lims = [min([scores_raw(mask); scores_loc(mask)]), max([scores_raw(mask); scores_loc(mask)])];
        plot(lims, lims, 'r--', 'LineWidth', 1);

        xlabel('scores\_raw (Sj)'); ylabel('scores\_loc (SjLoc)');
        title(sprintf('Sj vs SjLoc (order=%d, top %d%%)', ordr, 100-topP));
        axis equal; xlim([0 0.7]); ylim([0 0.8]); hold off;

        exportgraphics(figure(2), fullfile(outdir, sprintf('score_beforeafter_%s.png', tag)), 'Resolution', 300);
        savefig(2,              fullfile(outdir, sprintf('score_beforeafter_%s.fig', tag)));
      end

      % ===== Fig6: 시간 히스토그램(겹쳐, 선형축 그대로) =====
      has_time_raw = size(X,2) >= 3;
      has_time_loc = size(X,2) >= 4;
      if has_time_raw || has_time_loc
        raw_ms = []; loc_ms = [];
        if has_time_raw, raw_ms = X(:,3); raw_ms = raw_ms(isfinite(raw_ms)); end
        if has_time_loc, loc_ms = X(:,4); loc_ms = loc_ms(isfinite(loc_ms)); end

        figure(6); clf; hold on; grid on;
        if ~isempty(raw_ms)
          [~, e_r] = histcounts(raw_ms, 'BinMethod','fd');
          histogram(raw_ms, e_r, 'FaceAlpha',0.5,'EdgeAlpha',0.5);
        end
        if ~isempty(loc_ms)
          [~, e_l] = histcounts(loc_ms, 'BinMethod','fd');
          histogram(loc_ms, e_l, 'FaceAlpha',0.5,'EdgeAlpha',0.5);
        end
        xlabel('ms'); ylabel('count'); title(sprintf('Iteration time (order=%d)', ordr));
        xlim([0 timelim(oi)]); hold off;

        exportgraphics(figure(6), fullfile(outdir, sprintf('time_%s.png', tag)), 'Resolution', 300);
        savefig(6,              fullfile(outdir, sprintf('time_%s.fig', tag)));
      end

    end
  end
end
