PP01 = generateRandomPointsOnHexagonPrism(9800)+randn(9800,3)*0.007;
PP02 = generateRandomPointsOnSpaceship(4000)+randn(4000,3)*0.006;
PP03 = generateRandomPointsOnCylinder(8001)+randn(8001,3)*0.005;
PP04 = (2*randn(200,3)-1)*2;
PP05 = generateRandomPointsOnIcosahedron(4000);


PPm_use = [PP01 ; PP04];

%%
syms x y z 
order = 4;

TermsC = homogeneFischerTerms(order);
FuncsC = matlabFunction(TermsC);
Qinit = [1;0.0;0.0;0.00];
N = length(TermsC);

[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);
warning('off','MATLAB:nearlySingularMatrix')


%% 신형 ransac _ 속도 증가 by funcs 
%{
ransacPar = struct( ...
  'maxIter', 1000, ...
  'conf', 0.96, ...
  'thresh', 0.454, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0, ...   % 로컬 최적화 비활성화
  'locIters', 4, ... 
  'momentum', 0.66,...
  'damping', 0.85, ...
  'reg', 1e-6);
ransacPar.mc = struct( ...
    'on', true, ...
    'saveVarName', 'MC_scores_8', ...  % base workspace에 MC_scores 라는 변수로 저장
    'saveMatFile', '', ...
    'time', true, ...% 파일 저장 안 함
    'localOff', true);     % 로컬 완전 비활성화 (선택)
ransacPar.sim = struct( ...
    'freezeW',     true, ...  % w 갱신 건너뜀
    'freezeIter',  true, ...  % p.maxIter 동적 갱신 금지
    'freezeOmega', true ...   % omega 갱신 건너뜀
);
[DispRAN3, BetaRAN3, inlierMaskRAN3] = ...
  PoliNavigationSolver3_2_FischerRansac_MC(0, PPm_use, order, length(TermsC), Funcs1, Grads1, ransacPar);


if ~isempty(BetaRAN3)
    PPm_use_unbias = PPm_use - DispRAN2;
    
end
%% histogram
T=MC_scores_8;
scores_raw = T(:,1);
scores_loc = T(:,2);

raw = scores_raw(isfinite(scores_raw));
loc = scores_loc(isfinite(scores_loc));   % 로컬 미진입 row는 NaN이므로 제외

% 공통 bin 경계(두 분포 비교를 위해 동일한 bins 사용)
allvals = [raw; loc];
if isempty(allvals)
    error('MC_scores가 비어 있습니다.');
end
[~, edges] = histcounts(allvals, 'BinMethod','fd');  % Freedman–Diaconis
p95_raw = prctile(raw, 95);
p95_loc = NaN;
if ~isempty(loc)
    p95_loc = prctile(loc, 95);
end

figure(1)
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');
nexttile;
histogram(raw, edges); grid on; title('Sj (raw)');
xlabel('score'); ylabel('count');
hold on;
xline(p95_raw, 'r', 'LineWidth', 1.5);
yl = ylim; text(p95_raw, yl(2)*0.9, sprintf('95%% = %.4g', p95_raw), ...
    'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
hold off;

nexttile;
histogram(loc, edges); grid on; title('SjLoc (local)');
xlabel('score'); ylabel('count');
hold on;
xline(p95_loc, 'r', 'LineWidth', 1.5);
xl = xlim; yl = ylim;
xpos_loc = min(max(p95_loc, xl(1)), xl(2));               % 축 안으로 클램프
if xpos_loc > xl(1) + 0.9*(xl(2)-xl(1))                   % 오른쪽 가장자리면
    xpos_loc = xl(2) - 0.02*(xl(2)-xl(1));
    hal_loc = 'right';
else
    hal_loc = 'left';
end
text(xpos_loc, yl(2)*0.9, sprintf('95%% = %.4g', p95_loc), ... % ← p95_loc로 수정
    'Color','r','HorizontalAlignment',hal_loc,'VerticalAlignment','top');
hold off;

%시간 복잡도
raw_ms = T(:,3);
loc_ms = T(:,4);

raw_ms = raw_ms(isfinite(raw_ms));
loc_ms = loc_ms(isfinite(loc_ms));
[~, edges_raw] = histcounts(raw_ms, 'BinMethod','fd');
[~, edges_loc] = histcounts(loc_ms, 'BinMethod','fd');

figure(6); clf; hold on; grid on;
h1 = histogram(raw_ms, edges_raw, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
h2 = histogram(loc_ms, edges_loc, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
set(gca,'XScale','log');
xlabel('ms'); ylabel('count');
title('Iteration time overlay (raw\_ms vs loc\_ms)');
legend([h1 h2], {'raw\_ms','loc\_ms'}, 'Location','best');
hold off;
%%
% MC_scores_2: [Sj(raw), SjLoc(local)]
%scores_raw = MC_scores_4(:,1);
%scores_loc = MC_scores_4(:,2);

mask = isfinite(scores_raw) & isfinite(scores_loc);

% 상위 퍼센타일로 확대 기준 설정 (원하면 90, 97.5 등으로 조정)
topP = 90; topraw = 99.8; thresh = 0.6;
xr = prctile(scores_raw(mask), [topraw 100]);
yr = prctile(scores_loc(mask), [topP 100]);
xr = thresh * xr;
% 상위 구간에 속하는 점들만 강조
zoomMask = mask & scores_loc >= yr(1);
threshMask = mask & scores_raw >= xr(1);
figure(2); clf; hold on; grid on;

% 1) 전체 분포(연한 점)
plot(scores_raw(mask), scores_loc(mask), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6);

% 2) 상위 구간 강조(진한 점)
plot(scores_raw(zoomMask), scores_loc(zoomMask), '.', 'MarkerSize', 10);
plot(scores_raw(threshMask), scores_loc(threshMask), '.', 'MarkerSize', 5);

% 대각선 y=x
lims = [min([scores_raw(mask); scores_loc(mask)]), max([scores_raw(mask); scores_loc(mask)])];
plot(lims, lims, 'r--', 'LineWidth', 1);

% 축 범위: 상위 퍼센타일로 확대

%ylim([yr(1) yr(2)]);

xlabel('scores\_raw (Sj)');
ylabel('scores\_loc (SjLoc)');
title(sprintf('Sj vs SjLoc (zoomed to top %d%%)', 100-topP));
axis equal
hold off;
%}

%% 반복
totalN = 10000;
%eps_list      = [0.20, 0.30, 0.35, 0.40];
locIters_list = [2,4,7];
eps_list      = [0.00, 0.02, 0.05, 0.10];
%locIters_list = [2,4,7];
Niter         = 5000;
Thresh        = 0.300;
Thresh_list   = [0.2 0.3 0.45 0.6];
ransacPar = struct( ...
  'maxIter', Niter, ...
  'conf', 0.96, ...
  'thresh', Thresh, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0, ...   % 로컬 최적화 비활성화
  'locIters', 4, ... 
  'momentum', 0.66,...
  'damping', 0.85, ...
  'reg', 1e-6);
ransacPar.mc = struct( ...
    'on', true, ...
    'saveVarName', 'MC_scores', ...  % base workspace에 MC_scores 라는 변수로 저장
    'saveMatFile', '', ...
    'time', true, ...% 파일 저장 안 함
    'localOff', false);     % 로컬 완전 비활성화 (선택)
ransacPar.sim = struct( ...
    'freezeW',     true, ...  % w 갱신 건너뜀
    'freezeIter',  true, ...  % p.maxIter 동적 갱신 금지
    'freezeOmega', true ...   % omega 갱신 건너뜀
);
PP11 = generateRandomPointsOnHexagonPrism(10400)+randn(10400,3)*0.007;
PP14 = (2*randn(5000,3)-1)*2;

for e = eps_list
    n_out = round(totalN*e);  n_in = totalN - n_out;
    PPm_use = [ PP11(1:n_in,:); PP14(1:n_out,:) ];

    for L = locIters_list
        
            % ---- 여기서 '필드만' 수정 ----
            ransacPar.maxIter   = Niter;
            ransacPar.locIters  = L;
            ransacPar.mc.saveVarName = sprintf('MC_scores_w%02d_th%d_N%d_nT14', round(100*e), L, Niter);

            % 호출 (DisplacementLocal이 p를 받도록 내부에서 4번째 인자로 전달되어 있어야 함)
            [DispRAN, BetaRAN, inlierMaskRAN] = ...
              PoliNavigationSolver3_2_FischerRansac_MC(0, PPm_use, order, numel(TermsC), Funcs1, Grads1, ransacPar);
        
    end
end
%%
%BetaRAN = BetaRAN8;
%DispRAN = DispRAN8;
%inlierMaskRAN = inlierMaskRAN8;


if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN);

    PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = Funcs1_RAN(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,[-20 20 -10 10 -30 30],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 150);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-25.5 25.5 -10.5 10.5 -15.5 15.5])

    
    figure(3)
    scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,'b','filled');
    
    fimplicit3(f1_RAN-1,[-20 20 -20 20 -30 30],'FaceColor', [0.95, 0.82, 0.5],'EdgeColor','none','FaceAlpha',0.5, 'MeshDensity', 150);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-25.5 25.5 -15.5 15.5 -15.5 15.5])

end
%% histogram 바동
Niter = 5000;


% 대상 조합 (필요시 수정)
%w_list_pct = [0 2 5 10];     % w00, w02, w05, w10  → 4개
%w_list_pct = [20 30 35 40];     % w20, w30, w35, w40  → 4개
w_list_pct = [35 40];
L_list     = [2 4 7];        % L2, L4, L7          → 3개

% 저장 폴더
outdir = fullfile(pwd, 'image_ISS');
if ~exist(outdir,'dir'), mkdir(outdir); end

for wi = 1:numel(w_list_pct)
  for li = 1:numel(L_list)
    wpct = w_list_pct(wi);
    L    = L_list(li);

    % 변수명/태그 구성
    mcvar = sprintf('MC_scores_w%02d_L%d_N%d_nT14', wpct, L, Niter);
    tag   = erase(mcvar, 'MC_scores_');

    % 작업공간에 존재 확인 후 로드
    if ~evalin('base', sprintf('exist(''%s'',''var'')', mcvar))
        warning('없음: %s (건너뜀)', mcvar);
        continue;
    end
    X = evalin('base', mcvar);   % [Sj, SjLoc, raw_ms, loc_ms] (열 3,4는 없을 수 있음)

    % ----- Fig1: score 히스토그램 (공통 bins) -----



    scores_raw = X(:,1); m1 = isfinite(scores_raw);
    scores_loc = X(:,2); m2 = isfinite(scores_loc);

    if any(m1|m2)

       


        [~, edges] = histcounts([scores_raw(m1); scores_loc(m2)], 'BinMethod','fd');
        p95_score_raw = prctile(scores_raw(m1), 95);
        p95_score_loc = NaN;
        if ~isempty(loc)
            p95_score_loc = prctile(scores_loc(m2), 95);
        end
        figure(1)
        tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');
        nexttile;
        histogram(scores_raw, edges); grid on; title('Sj (raw)');
        xlabel('score'); ylabel('count');
        hold on;
        xline(p95_score_raw, 'r', 'LineWidth', 1.5);
        yl = ylim; text(p95_score_raw, yl(2)*0.9, sprintf('95%% = %.4g', p95_raw), ...
            'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
        hold off;
        
        nexttile;
        histogram(scores_loc, edges); grid on; title('SjLoc (local)');
        xlabel('score'); ylabel('count');
        hold on;
        xline(p95_score_loc, 'r', 'LineWidth', 1.5);
        xl = xlim; yl = ylim;
        xpos_loc = min(max(p95_score_loc, xl(1)), xl(2));               % 축 안으로 클램프
        if xpos_loc > xl(1) + 0.9*(xl(2)-xl(1))                   % 오른쪽 가장자리면
            xpos_loc = xl(2) - 0.02*(xl(2)-xl(1));
            hal_loc = 'right';
        else
            hal_loc = 'left';
        end
        text(xpos_loc, yl(2)*0.9, sprintf('95%% = %.4g', p95_score_loc), ... % ← p95_loc로 수정
            'Color','r','HorizontalAlignment',hal_loc,'VerticalAlignment','top');
        hold off;
        exportgraphics(figure(1), fullfile(outdir, sprintf('score_%s.png', tag)), 'Resolution', 300);
        savefig(1,              fullfile(outdir, sprintf('score_%s.fig', tag)));
    end

    % ----- Fig2: Sj vs SjLoc 산점도 -----
    mask = isfinite(scores_raw) & isfinite(scores_loc);
    topP = 90; topraw = 99.8; thresh = 0.6;
    xr = prctile(scores_raw(mask), [topraw 100]);
    yr = prctile(scores_loc(mask), [topP 100]);
    xr = thresh * xr;
    % 상위 구간에 속하는 점들만 강조
    zoomMask = mask & scores_loc >= yr(1);
    threshMask = mask & scores_raw >= xr(1);
    figure(2); clf; hold on; grid on;
    
    % 1) 전체 분포(연한 점)
    plot(scores_raw(mask), scores_loc(mask), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6);
    
    % 2) 상위 구간 강조(진한 점)
    plot(scores_raw(zoomMask), scores_loc(zoomMask), '.', 'MarkerSize', 10);
    plot(scores_raw(threshMask), scores_loc(threshMask), '.', 'MarkerSize', 5);
    
    % 대각선 y=x
    lims = [min([scores_raw(mask); scores_loc(mask)]), max([scores_raw(mask); scores_loc(mask)])];
    plot(lims, lims, 'r--', 'LineWidth', 1);
    
    % 축 범위: 상위 퍼센타일로 확대
    
    %ylim([yr(1) yr(2)]);
    outdir = fullfile(pwd, 'image_ISS');
    xlabel('scores\_raw (Sj)');
    ylabel('scores\_loc (SjLoc)');
    title(sprintf('Sj vs SjLoc (zoomed to top %d%%)', 100-topP));
    axis equal
    hold off;
    exportgraphics(figure(2), fullfile(outdir, sprintf('score_beforeafter_%s.png', tag)), 'Resolution', 300);
    savefig(2,              fullfile(outdir, sprintf('score_beforeafter_%s.fig', tag)));

    % ----- Fig6: 시간 히스토그램(겹쳐) -----
    has_time_raw = size(X,2) >= 3;
    has_time_loc = size(X,2) >= 4;
    if has_time_raw || has_time_loc
        raw_ms = []; loc_ms = [];
        if has_time_raw, raw_ms = X(:,3); raw_ms = raw_ms(isfinite(raw_ms)); end
        if has_time_loc, loc_ms = X(:,4); loc_ms = loc_ms(isfinite(loc_ms)); end

        figure(6); clf; hold on; grid on;
        if ~isempty(raw_ms)
            [~, e_r] = histcounts(raw_ms, 'BinMethod','fd');
            h1 = histogram(raw_ms, e_r, 'FaceAlpha',0.5,'EdgeAlpha',0.5);
        end
        if ~isempty(loc_ms)
            [~, e_l] = histcounts(loc_ms, 'BinMethod','fd');
            h2 = histogram(loc_ms, e_l, 'FaceAlpha',0.5,'EdgeAlpha',0.5);
        end
        xlabel('ms'); ylabel('count'); title('Iteration time overlay');
        lg = {};
        if exist('h1','var'), lg{end+1}='raw\_ms'; end
        if exist('h2','var'), lg{end+1}='loc\_ms'; end
        if ~isempty(lg), legend(lg,'Location','best'); end
        % x축 범위 양쪽 데이터 커버
        xl = xlim;
        if exist('e_r','var'), xl(1)=min(xl(1), e_r(1));  xl(2)=max(xl(2), e_r(end)); end
        if exist('e_l','var'), xl(1)=min(xl(1), e_l(1));  xl(2)=max(xl(2), e_l(end)); end
        xlim(xl);
        hold off;

        exportgraphics(figure(6), fullfile(outdir, sprintf('time_%s.png', tag)), 'Resolution', 300);
        savefig(6,              fullfile(outdir, sprintf('time_%s.fig', tag)));
    end
  end
end

%%
%{
mcvar = 'MC_scores_w00_L2_N1000_nT10';
%T=MC_scores_w00_L2_N1000;

if exist('mcvar','var') && ~isempty(mcvar)
    if evalin('base', sprintf('exist(''%s'',''var'')', mcvar))
        T   = evalin('base', mcvar);           % 데이터 읽기
        tag = erase(mcvar, 'MC_scores_');      % 저장용 접미사
    else
        error('변수 %s 가 workspace에 없습니다.', mcvar);
    end
elseif exist('T','var') && ~isempty(T)
    
    if ~exist('tag','var') || isempty(tag), tag = 'custom'; end
end

scores_raw = T(:,1);
scores_loc = T(:,2);

raw = scores_raw(isfinite(scores_raw));
loc = scores_loc(isfinite(scores_loc));   % 로컬 미진입 row는 NaN이므로 제외

% 공통 bin 경계(두 분포 비교를 위해 동일한 bins 사용)
allvals = [raw; loc];
if isempty(allvals)
    error('MC_scores가 비어 있습니다.');
end
[~, edges] = histcounts(allvals, 'BinMethod','fd');  % Freedman–Diaconis
p95_raw = prctile(raw, 95);
p95_loc = NaN;
if ~isempty(loc)
    p95_loc = prctile(loc, 95);
end

figure(1)
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');
nexttile;
histogram(raw, edges); grid on; title('Sj (raw)');
xlabel('score'); ylabel('count');
hold on;
xline(p95_raw, 'r', 'LineWidth', 1.5);
yl = ylim; text(p95_raw, yl(2)*0.9, sprintf('95%% = %.4g', p95_raw), ...
    'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
hold off;

nexttile;
histogram(loc, edges); grid on; title('SjLoc (local)');
xlabel('score'); ylabel('count');
hold on;
xline(p95_loc, 'r', 'LineWidth', 1.5);
xl = xlim; yl = ylim;
xpos_loc = min(max(p95_loc, xl(1)), xl(2));               % 축 안으로 클램프
if xpos_loc > xl(1) + 0.9*(xl(2)-xl(1))                   % 오른쪽 가장자리면
    xpos_loc = xl(2) - 0.02*(xl(2)-xl(1));
    hal_loc = 'right';
else
    hal_loc = 'left';
end
text(xpos_loc, yl(2)*0.9, sprintf('95%% = %.4g', p95_loc), ... % ← p95_loc로 수정
    'Color','r','HorizontalAlignment',hal_loc,'VerticalAlignment','top');
hold off;

%시간 복잡도
raw_ms = T(:,3);
loc_ms = T(:,4);

raw_ms = raw_ms(isfinite(raw_ms));
loc_ms = loc_ms(isfinite(loc_ms));
[~, edges_raw] = histcounts(raw_ms, 'BinMethod','fd');
[~, edges_loc] = histcounts(loc_ms, 'BinMethod','fd');

figure(6); clf; hold on; grid on;
h1 = histogram(raw_ms, edges_raw, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
h2 = histogram(loc_ms, edges_loc, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
set(gca,'XScale','log');
xlabel('ms'); ylabel('count');
title('Iteration time overlay (raw\_ms vs loc\_ms)');
legend([h1 h2], {'raw\_ms','loc\_ms'}, 'Location','best');
hold off;

mask = isfinite(scores_raw) & isfinite(scores_loc);
topP = 90; topraw = 99.8; thresh = 0.6;
xr = prctile(scores_raw(mask), [topraw 100]);
yr = prctile(scores_loc(mask), [topP 100]);
xr = thresh * xr;
% 상위 구간에 속하는 점들만 강조
zoomMask = mask & scores_loc >= yr(1);
threshMask = mask & scores_raw >= xr(1);
figure(2); clf; hold on; grid on;

% 1) 전체 분포(연한 점)
plot(scores_raw(mask), scores_loc(mask), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 6);

% 2) 상위 구간 강조(진한 점)
plot(scores_raw(zoomMask), scores_loc(zoomMask), '.', 'MarkerSize', 10);
plot(scores_raw(threshMask), scores_loc(threshMask), '.', 'MarkerSize', 5);

% 대각선 y=x
lims = [min([scores_raw(mask); scores_loc(mask)]), max([scores_raw(mask); scores_loc(mask)])];
plot(lims, lims, 'r--', 'LineWidth', 1);

% 축 범위: 상위 퍼센타일로 확대

%ylim([yr(1) yr(2)]);
outdir = fullfile(pwd, 'image_ISS');
xlabel('scores\_raw (Sj)');
ylabel('scores\_loc (SjLoc)');
title(sprintf('Sj vs SjLoc (zoomed to top %d%%)', 100-topP));
axis equal
hold off;
exportgraphics(figure(1), fullfile(outdir, sprintf('score_%s.png',             tag)), 'Resolution', 300);
savefig(1,              fullfile(outdir, sprintf('score_%s.fig',             tag)));

exportgraphics(figure(2), fullfile(outdir, sprintf('score_beforeafter_%s.png', tag)), 'Resolution', 300);
savefig(2,              fullfile(outdir, sprintf('score_beforeafter_%s.fig', tag)));

exportgraphics(figure(6), fullfile(outdir, sprintf('time_%s.png',              tag)), 'Resolution', 300);
savefig(6,              fullfile(outdir, sprintf('time_%s.fig',
tag)));
%}

