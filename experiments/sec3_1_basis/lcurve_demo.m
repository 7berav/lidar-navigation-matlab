%% lcurve_demo.m — L-curve + CD/Hausdorff(occluded) 검증 (Sec 3.1b)
%
% 세팅:
%   occlusion : 도형을 뷰마다 무작위 회전 → z 상위 keepFrac(=1-occ)만 visible
%   스윕      : occ ∈ {0.30,0.50,0.70}, nViews=10 (뷰 기하 통계), Sobolev H^1 ridge
%   피팅      : visible에서 RANSAC 드로우 크기 N_fit = k×28 (k=1.0/1.2/1.4) 무작위 서브샘플
%   L-curve   : 서브샘플 대수 MSE (y) vs 페널티 seminorm β'Dβ (x)
%   검증      : hidden GT 점 ↔ 등위면 f=1 의 Chamfer(CD_occ) + Hausdorff95(HD_occ)
%   메시지    : ridge 이득이 뷰 기하에 강건하고 occlusion↑ 일수록 커짐
%
% Fischer 6차 28항, 이상치 없음(eps=0). 실행: MATLAB R2025a (Symbolic+Statistics).

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
rng(7);
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

% 배치 실행 전 workspace에서 덮어쓸 수 있는 실험 설정.
% 예: shapeName='cube'; quickMode=true; run('.../lcurve_demo.m')
if ~exist('shapeName', 'var') || isempty(shapeName), shapeName = 'hexagon'; end
if ~exist('quickMode', 'var') || isempty(quickMode), quickMode = false; end
if ~exist('kFocus', 'var') || isempty(kFocus), kFocus = 1.4; end
if ~exist('order', 'var') || isempty(order), order = 4; end
shapeName = lower(char(shapeName));

% ---- 출력 폴더: 기존 규약대로 날짜 스탬프 ----
outDir   = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'test_ridge');
if ~isfolder(outDir), mkdir(outDir); end
modeTag  = '';
if quickMode, modeTag = '_quick'; end
runStamp = sprintf('%s_%s_or%d%s', datestr(now, 'yyyymmdd_HHMM'), shapeName, order, modeTag);
figDir   = fullfile(outDir, runStamp);
mkdir(figDir);
fprintf('그림 폴더 → %s  (order=%d)\n', figDir, order);

%% ---- 기저 ----
TermsC = homogeneFischerTerms(order);
[Funcs1, ~, nT] = makeFuncsGradsStack(TermsC);

% Sobolev 단위 페널티 가중 (regressionFourthOrder 규약과 동일)
levels = order:-2:0;
wSob   = zeros(nT, 1);
ix = 1;
for l = levels
    nk = 2*l + 1;
    wSob(ix:ix+nk-1) = (l*(l+1)) / (order*(order+1));
    ix = ix + nk;
end

%% ---- 등위면 샘플용 그리드 ----
gLin = linspace(-1.6, 1.6, 49);
[GX, GY, GZ] = meshgrid(gLin, gLin, gLin);
PhiGrid = calculateFourthOrder([GX(:), GY(:), GZ(:)], Funcs1);

%% ---- 스윕 파라미터 ----
totalN    = 6000;
noiseSig  = 0.007;
if quickMode
    occ_list = [0.30, 0.50, 0.70];
    nViews   = 3;
    k_list   = kFocus;
    nDrawScal= 3;
else
    occ_list = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70];
    nViews   = 10;
    k_list   = [1.0, 1.2, 1.4];
    nDrawScal= 10;
end
nOcc      = numel(occ_list);
snapOcc   = intersect([0.30, 0.50, 0.70], occ_list, 'stable');
panelOcc  = intersect([0.30, 0.50, 0.70], occ_list, 'stable');
nK        = numel(k_list);
% CD_all 최소가 3e-1~1e0 부근에 있어(광역 스윕 확인) 1e0, 3e0까지 감싼다.
% 3e0 초과는 df<3 구-붕괴 점근구간이라 제외.
lambdas   = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0];
nL        = numel(lambdas);
hdPctl    = 95;                      % 로버스트 Hausdorff 백분위수
nSurfEval = 15000;                   % mesh 삼각형 면적비례 평가점 수

[~, kFocusIdx] = min(abs(k_list - kFocus));
kFocus = k_list(kFocusIdx);          % quick/full 모두 실제 k_list 원소에 맞춤

% 베이스 점군 1회 생성. 같은 10개 회전/노이즈 점군을 모든 occlusion에서 공유한다.
Pbase = makeShapePoints(shapeName, totalN);
PcleanViews = zeros(totalN, 3, nViews);
PnoisyViews = zeros(totalN, 3, nViews);
for iv = 1:nViews
    axisv = randn(3,1); axisv = axisv / norm(axisv);
    ang   = 2*pi*rand;
    K = [0 -axisv(3) axisv(2); axisv(3) 0 -axisv(1); -axisv(2) axisv(1) 0];
    Rv = eye(3) + sin(ang)*K + (1-cos(ang))*(K*K);
    PcleanViews(:,:,iv) = Pbase * Rv.';
    PnoisyViews(:,:,iv) = PcleanViews(:,:,iv) + noiseSig * randn(totalN, 3);
end

% 결과 저장: [occ, view, k, lambda]
MSE = nan(nOcc, nViews, nK, nL);  PEN = nan(nOcc, nViews, nK, nL);
GCVv= nan(nOcc, nViews, nK, nL);  CDO = nan(nOcc, nViews, nK, nL);
HDO = nan(nOcc, nViews, nK, nL);  CDA = nan(nOcc, nViews, nK, nL);
HDA = nan(nOcc, nViews, nK, nL);  DFm = nan(nOcc, nViews, nK, nL);
OSF = nan(nOcc, nViews, nK, nL);  % sampled fitted surface 중 hidden 영역 비율
DFraw = nan(nOcc, nViews, nK, nL, nDrawScal);

% 전체 visible baseline: [occ, view, lambda]
MSEF = nan(nOcc, nViews, nL);  PENF = nan(nOcc, nViews, nL);
GCVF = nan(nOcc, nViews, nL);  DFF  = nan(nOcc, nViews, nL);
CDOF = nan(nOcc, nViews, nL);  HDOF = nan(nOcc, nViews, nL);
CDAF = nan(nOcc, nViews, nL);  HDAF = nan(nOcc, nViews, nL);
OSFF = nan(nOcc, nViews, nL);
% fig3 스냅샷용: occ=0.50/0.70, view 1, 모든 k와 lambda의 β
snapBeta = cell(nOcc, nK, nL);
snapBetaFull = cell(nOcc, nL);
snapVis  = cell(nOcc, 1);
snapHid  = cell(nOcc, 1);
snapCen  = cell(nOcc, nK);   % 스냅 드로우 무게중심(등위면 원프레임 복원용)
snapCenFull = cell(nOcc, 1);  % 전체 visible baseline 중심
snapFit  = cell(nOcc, nK);   % fig3 모델이 실제로 사용한 rep1 부분집합

Trows = {};   % long-form CSV 누적

for io = 1:nOcc
    occ      = occ_list(io);
    keepFrac = 1 - occ;
    for iv = 1:nViews
        % 동일한 view 번호는 모든 occlusion에서 같은 회전과 노이즈를 사용
        Pclean = PcleanViews(:,:,iv);
        Pnoisy = PnoisyViews(:,:,iv);
        zThr   = prctile(Pclean(:,3), 100*(1 - keepFrac));
        vis    = Pclean(:,3) >= zThr;
        Pvis   = Pnoisy(vis, :);
        PallClean = Pclean;
        PhidClean = Pclean(~vis, :);
        emptyPenalty = norm(max(PallClean,[],1)-min(PallClean,[],1));

        NvisAll = size(Pvis, 1);

        % ---- 전체 visible baseline (hidden GT는 피팅에 절대 사용하지 않음) ----
        cFull   = mean(Pvis, 1);
        PhiFull = calculateFourthOrder(Pvis - cFull, Funcs1);
        A0Full  = PhiFull.' * PhiFull;
        rhsFull = PhiFull.' * ones(NvisAll, 1);
        for i = 1:nL
            lam   = lambdas(i);
            DFull = lam * NvisAll * diag(wSob);
            AFull = A0Full + DFull;
            bFull = AFull \ rhsFull;
            rFull = PhiFull*bFull - 1;

            MSEF(io,iv,i) = mean(rFull.^2);
            PENF(io,iv,i) = bFull.'*(wSob.*bFull);
            DFF(io,iv,i)  = trace(AFull \ A0Full);
            GCVF(io,iv,i) = MSEF(io,iv,i) / max(1 - DFF(io,iv,i)/NvisAll, 1e-12)^2;

            Vall = meshSurfacePoints(bFull, PhiGrid, GX, GY, GZ, cFull, nSurfEval);
            Vocc = Vall(Vall(:,3) < zThr, :);
            OSFF(io,iv,i) = size(Vocc,1)/max(size(Vall,1),1);
            [CDAF(io,iv,i), HDAF(io,iv,i)] = symmetricMetrics(PallClean,Vall,hdPctl,emptyPenalty);
            [CDOF(io,iv,i), HDOF(io,iv,i)] = symmetricMetrics(PhidClean,Vocc,hdPctl,emptyPenalty);

            if any(abs(occ-snapOcc)<1e-9) && iv==1
                snapBetaFull{io,i} = bFull;
                snapCenFull{io} = cFull;
            end

            Trows(end+1,:) = {shapeName, 'full', occ, iv, NaN, NvisAll, lam, ...
                MSEF(io,iv,i), PENF(io,iv,i), GCVF(io,iv,i), DFF(io,iv,i), ...
                CDOF(io,iv,i), HDOF(io,iv,i), CDAF(io,iv,i), HDAF(io,iv,i),OSFF(io,iv,i)}; %#ok<AGROW>
        end

        for kk = 1:nK
            Nfit = round(k_list(kk) * nT);
            % 스칼라용 서브샘플 풀 (동일 드로우를 λ 간 공유)
            poolIdx = zeros(nDrawScal, Nfit);
            for rr = 1:nDrawScal
                poolIdx(rr, :) = randperm(NvisAll, Nfit);
            end
            % 드로우별 무게중심 재센터링(솔버 :166 규약). Φ는 λ 무관 → λ 루프 밖에서 캐시.
            cPool   = zeros(nDrawScal, 3);
            PhiPool = cell(nDrawScal, 1);
            for rr = 1:nDrawScal
                cPool(rr,:) = mean(Pvis(poolIdx(rr,:), :), 1);
                PhiPool{rr} = calculateFourthOrder( ...
                    Pvis(poolIdx(rr,:), :) - cPool(rr,:), Funcs1);
            end
            % 거리 지표는 rep1 드로우 사용 → PhiPool{1}, cPool(1,:)
            cCD   = cPool(1, :);            % rep1 드로우 무게중심(등위면 원프레임 복원용)

            for i = 1:nL
                lam = lambdas(i);
                D   = lam * Nfit * diag(wSob);

                % --- 스칼라 지표: nDrawScal 중앙값 ---
                mseR = zeros(nDrawScal,1); penR = mseR; gcvR = mseR; dfR = mseR;
                for rr = 1:nDrawScal
                    Phi_s = PhiPool{rr};
                    A0 = Phi_s.'*Phi_s;  A = A0 + D;
                    b  = A \ (Phi_s.'*ones(Nfit,1));
                    r  = Phi_s*b - 1;
                    mseR(rr) = mean(r.^2);
                    penR(rr) = b.'*(wSob.*b);
                    trH = trace(A \ A0);
                    dfR(rr)  = trH;                    % df(λ)=tr(H_λ)
                    gcvR(rr) = mseR(rr) / max(1 - trH/Nfit, 1e-12)^2;
                end
                MSE(io,iv,kk,i) = median(mseR);
                PEN(io,iv,kk,i) = median(penR);
                GCVv(io,iv,kk,i)= median(gcvR);
                DFm(io,iv,kk,i) = median(dfR);
                DFraw(io,iv,kk,i,:) = dfR;

                % --- 거리 지표: rep1 mesh를 면적비례 샘플링한 clean-GT all/occ 평가 ---
                Phi_s = PhiPool{1};
                bCD = (Phi_s.'*Phi_s + D) \ (Phi_s.'*ones(Nfit,1));
                Vall = meshSurfacePoints(bCD, PhiGrid, GX, GY, GZ, cCD, nSurfEval);
                Vocc = Vall(Vall(:,3) < zThr, :);
                OSF(io,iv,kk,i) = size(Vocc,1)/max(size(Vall,1),1);
                [CDA(io,iv,kk,i), HDA(io,iv,kk,i)] = symmetricMetrics(PallClean,Vall,hdPctl,emptyPenalty);
                [CDO(io,iv,kk,i), HDO(io,iv,kk,i)] = symmetricMetrics(PhidClean,Vocc,hdPctl,emptyPenalty);

                Trows(end+1,:) = {shapeName, 'subset', occ, iv, k_list(kk), Nfit, lam, ...
                    MSE(io,iv,kk,i), PEN(io,iv,kk,i), GCVv(io,iv,kk,i), ...
                    DFm(io,iv,kk,i), CDO(io,iv,kk,i), HDO(io,iv,kk,i), ...
                    CDA(io,iv,kk,i), HDA(io,iv,kk,i),OSF(io,iv,kk,i)}; %#ok<AGROW>

                if any(abs(occ - snapOcc) < 1e-9) && iv == 1
                    snapBeta{io,kk,i} = bCD;
                    snapCen{io,kk}    = cCD;
                    snapFit{io,kk}    = Pvis(poolIdx(1,:),:);
                    if isempty(snapVis{io})
                        snapVis{io} = Pvis;
                        snapHid{io} = PhidClean;
                    end
                end
            end
        end
        fprintf('occ=%.2f view=%2d/%d  (vis %d / hid %d)\n', ...
            occ, iv, nViews, size(Pvis,1), size(PhidClean,1));
    end
end

% 뷰 집계
mCDO = reshape(mean(CDO,2,'omitnan'), [nOcc,nK,nL]); sCDO = reshape(std(CDO,0,2,'omitnan'), [nOcc,nK,nL]);
mHDO = reshape(mean(HDO,2,'omitnan'), [nOcc,nK,nL]); sHDO = reshape(std(HDO,0,2,'omitnan'), [nOcc,nK,nL]);
mCDA = reshape(mean(CDA,2,'omitnan'), [nOcc,nK,nL]); sCDA = reshape(std(CDA,0,2,'omitnan'), [nOcc,nK,nL]);
mHDA = reshape(mean(HDA,2,'omitnan'), [nOcc,nK,nL]); sHDA = reshape(std(HDA,0,2,'omitnan'), [nOcc,nK,nL]);
mMSE = reshape(mean(MSE,2,'omitnan'), [nOcc,nK,nL]);
mPEN = reshape(mean(PEN,2,'omitnan'), [nOcc,nK,nL]);
mGCV = reshape(mean(GCVv,2,'omitnan'),[nOcc,nK,nL]);
mDF  = reshape(mean(DFm,2,'omitnan'), [nOcc,nK,nL]);
sDF  = reshape(std(DFm,0,2,'omitnan'),[nOcc,nK,nL]); % 뷰별 draw 중앙값의 표준편차

mMSEF = reshape(mean(MSEF,2,'omitnan'),[nOcc,nL]); mPENF = reshape(mean(PENF,2,'omitnan'),[nOcc,nL]);
mDFF  = reshape(mean(DFF,2,'omitnan'), [nOcc,nL]); sDFF  = reshape(std(DFF,0,2,'omitnan'),[nOcc,nL]);
mCDOF = reshape(mean(CDOF,2,'omitnan'),[nOcc,nL]); sCDOF = reshape(std(CDOF,0,2,'omitnan'),[nOcc,nL]);
mHDOF = reshape(mean(HDOF,2,'omitnan'),[nOcc,nL]); sHDOF = reshape(std(HDOF,0,2,'omitnan'),[nOcc,nL]);
mCDAF = reshape(mean(CDAF,2,'omitnan'),[nOcc,nL]); sCDAF = reshape(std(CDAF,0,2,'omitnan'),[nOcc,nL]);
mHDAF = reshape(mean(HDAF,2,'omitnan'),[nOcc,nL]); sHDAF = reshape(std(HDAF,0,2,'omitnan'),[nOcc,nL]);

%% ---- CSV ----
T = cell2table(Trows, 'VariableNames', ...
    {'shape','fit_type','occ','view','k','n_fit','lambda','mse','pen','gcv','df', ...
     'cd_occ','hd_occ','cd_all','hd_all','surface_occ_fraction'});
writetable(T, fullfile(figDir, 'lcurve_demo_results.csv'));
save(fullfile(figDir,'lcurve_demo_df_raw.mat'), ...
    'DFraw','DFm','DFF','occ_list','k_list','lambdas','shapeName', ...
    'snapBeta','snapBetaFull','snapCen','snapCenFull','snapFit','snapVis','snapHid');

%% ---- 그림 축 준비 ----
lamPlot = lambdas; lamPlot(1) = 3e-5;   % λ=0 로그축 표시용
colsK   = lines(nK);
panelIdx = find(ismember(round(occ_list*100), round(panelOcc*100)));  % fig1/2/5용 occ 인덱스
nPanel   = numel(panelIdx);

%% ---- Figure 1: L-curve (occ별 subplot, 모든 k, 뷰 평균+스프레드) ----
% k=1.0(N=nT=28)은 λ=0에서 보간→MSE≈0 이상치라 λ=0 점을 제외하고 그린다.
% y축 하한은 1e-2로 고정(의미 있는 점 min≈2.8e-2은 모두 유지, 보간 이상치만 배제).
yTop = max([mMSE(:); mMSEF(:)]) * 1.5;
fig = figure('Visible','off','Position',[50 50 1400 460]);
for p = 1:nPanel
    io = panelIdx(p);
    subplot(1, nPanel, p); hold on; grid on;
    for kk = 1:nK
        % k=1.0은 λ=0(보간 이상치) 제외
        if k_list(kk) == 1.0, useI = 2:nL; else, useI = 1:nL; end
        % 같은 k의 view별 곡선은 옅게, view 평균은 진하게 표시
        for iv = 1:nViews
            loglog(squeeze(PEN(io,iv,kk,useI)), squeeze(MSE(io,iv,kk,useI)), '-', ...
                'Color', [colsK(kk,:) 0.18], 'HandleVisibility','off');
        end
        penM = squeeze(mPEN(io,kk,useI)); mseM = squeeze(mMSE(io,kk,useI));
        jElb = mengerElbow(log(penM.'), log(mseM.'));
        iElb = useI(jElb);   % 전체 lambdas 인덱스로 복원
        loglog(penM, mseM, 'o-', 'Color', colsK(kk,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('k=%.1f (elbow %.0e)', k_list(kk), lambdas(iElb)));
        loglog(mPEN(io,kk,iElb), mMSE(io,kk,iElb), 's', 'Color', colsK(kk,:), ...
            'MarkerFaceColor', colsK(kk,:), 'MarkerSize', 9, ...
            'LineWidth', 1.5, 'HandleVisibility','off');
        if kk == kFocusIdx
            for i = [2 4 6 8]
                text(mPEN(io,kk,i)*1.06, mMSE(io,kk,i), sprintf('%.0e',lambdas(i)), ...
                    'Color', colsK(kk,:), 'FontSize', 7);
            end
        end
    end
    % 모든 visible 점으로 적합한 대조군
    penF = mPENF(io,:); mseF = mMSEF(io,:);
    iElbF = mengerElbow(log(penF), log(mseF));
    loglog(penF, mseF, 'k--', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('all visible (elbow %.0e)', lambdas(iElbF)));
    loglog(penF(iElbF), mseF(iElbF), 'ks', 'MarkerFaceColor','k', ...
        'MarkerSize', 9, 'HandleVisibility','off');
    set(gca,'XScale','log','YScale','log');
    ylim([1e-2, yTop]);
    xlabel('\beta^T W \beta'); ylabel('algebraic MSE');
    title(sprintf('occ=%d%% (vis %d%%)', ...
        round(100*occ_list(io)), round(100*(1-occ_list(io)))));
    legend('Location','best','FontSize',7);
end
sgtitle(sprintf('%s L-curve — color: subset k, dashed black: all visible, %d views', ...
    shapeName, nViews));
saveFigBoth(fig, figDir, 'fig1_lcurve');

%% ---- Figure 2: clean-GT occ/all 거리 지표 (각각 2×nPanel) ----
scopeNames = {'occ','all'};
metricNames = {{'CD_{occ}','HD_{occ} (95pct)'}, {'CD_{all}','HD_{all} (95pct)'}};
mSets = {{mCDO,mHDO}, {mCDA,mHDA}};  sSets = {{sCDO,sHDO}, {sCDA,sHDA}};
mFullSets = {{mCDOF,mHDOF}, {mCDAF,mHDAF}};
sFullSets = {{sCDOF,sHDOF}, {sCDAF,sHDAF}};
for sc = 1:2
    fig = figure('Visible','off','Position',[40 40 1400 760]);
    for mrow = 1:2
        for p = 1:nPanel
            io = panelIdx(p);
            subplot(2, nPanel, (mrow-1)*nPanel + p); hold on; grid on;
            for kk = 1:nK
                mu = squeeze(mSets{sc}{mrow}(io,kk,:));
                sd = squeeze(sSets{sc}{mrow}(io,kk,:));
                mu = mu(:); sd = sd(:);
                fill([lamPlot fliplr(lamPlot)], ...
                    [max(mu-sd,0).' fliplr((mu+sd).')], colsK(kk,:), ...
                    'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
                plot(lamPlot, mu, 'o-', 'Color',colsK(kk,:), 'LineWidth',1.4, ...
                    'DisplayName',sprintf('k=%.1f',k_list(kk)));
                [~,imin] = min(mu);
                plot(lamPlot(imin),mu(imin),'p','Color',colsK(kk,:), ...
                    'MarkerFaceColor',colsK(kk,:),'MarkerSize',9,'HandleVisibility','off');
            end
            muF = squeeze(mFullSets{sc}{mrow}(io,:));
            sdF = squeeze(sFullSets{sc}{mrow}(io,:));
            muF = muF(:); sdF = sdF(:);
            fill([lamPlot fliplr(lamPlot)], ...
                [max(muF-sdF,0).' fliplr((muF+sdF).')], [0.35 0.35 0.35], ...
                'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
            plot(lamPlot,muF,'k--','LineWidth',1.8,'DisplayName','all visible fit');
            [~,iminF] = min(muF);
            plot(lamPlot(iminF),muF(iminF),'kp','MarkerFaceColor','k', ...
                'MarkerSize',9,'HandleVisibility','off');

            set(gca,'XScale','log'); yl = ylim; ylim([0 yl(2)]);
            if mrow==2, xlabel('\lambda (OLS at 3e-5)'); end
            ylabel(metricNames{sc}{mrow});
            if mrow==1, title(sprintf('occ=%d%%',round(100*occ_list(io)))); end
            if p==1 && mrow==1, legend('Location','best','FontSize',7); end
        end
    end
    sgtitle(sprintf('%s %s-region clean-GT error — mean\\pmstd over %d views', ...
        shapeName, scopeNames{sc}, nViews));
    saveFigBoth(fig, figDir, sprintf('fig2_cd_hd_%s',scopeNames{sc}));
end

%% ---- Figure 3: 등위면 스냅샷 (snapOcc=30/50/70%, 모든 k, view1) ----
% 각 그림은 OLS, lambda=1e-3, lambda=1e0(CD-최적 부근)의 3-panel 비교다.
iShow = [find(lambdas==0), find(lambdas==1e-2), find(lambdas==1e0)];
for io = 1:nOcc
    if ~any(abs(occ_list(io) - snapOcc) < 1e-9), continue; end   % 대표 occ만
    for kk = 1:nK
        fig = figure('Visible','off','Position',[50 50 1400 460]);
        for s = 1:numel(iShow)
            i = iShow(s);
            subplot(1, numel(iShow), s); hold on;
            F = reshape(PhiGrid*snapBeta{io,kk,i}, size(GX));
            fvS = isosurface(GX, GY, GZ, F, 1);
            fvS.vertices = fvS.vertices + snapCen{io,kk};   % 센터링→원프레임 복원
            pa = patch(fvS);
            set(pa,'FaceColor',[0.2 0.6 0.9],'EdgeColor','none','FaceAlpha',0.45);
            Pvis = snapVis{io}; Phid = snapHid{io};
            cPlot = snapCen{io,kk};
            plot3(Pvis(1:4:end,1),Pvis(1:4:end,2),Pvis(1:4:end,3),'.', ...
                'Color',[0 0.3 0.8],'MarkerSize',3);
            plot3(Phid(1:4:end,1),Phid(1:4:end,2),Phid(1:4:end,3),'.', ...
                'Color',[0.6 0.6 0.6],'MarkerSize',3);
            plot3(cPlot(1),cPlot(2),cPlot(3),'mx', ...
                'MarkerSize',11,'LineWidth',2);
            axis equal; grid on; view(135,20); camlight; lighting gouraud;
            xlim([-1.6 1.6]); ylim([-1.6 1.6]); zlim([-1.6 1.6]);
            if lambdas(i)==0
                tt='OLS (\lambda=0)';
            else
                tt=sprintf('\\lambda=%.0e',lambdas(i));
            end
            % 표시하는 view1 표면과 동일한 CD를 제목에 사용한다.
            title(sprintf('%s  CD_{occ}=%.3f', tt, CDO(io,1,kk,i)));
        end
        Nfit = round(k_list(kk) * nT);
        sgtitle(sprintf(['occ=%d%%, k=%.1f (N=%d), view 1 — ' ...
            'blue: visible, gray: hidden clean GT, magenta: fit center'], ...
            round(100*occ_list(io)), k_list(kk), Nfit));
        figName = sprintf('fig3_o%02d_k%02d', ...
            round(100*occ_list(io)), round(10*k_list(kk)));
        saveFigBoth(fig, figDir, figName);
        close(fig);
    end
end

%% ---- Figure 3 full: 전체 visible baseline (동일 occ/view/lambda) ----
for io = 1:nOcc
    if ~any(abs(occ_list(io)-snapOcc)<1e-9), continue; end
    Pvis = snapVis{io}; Phid = snapHid{io}; cFull = snapCenFull{io};
    fig = figure('Visible','off','Position',[50 50 1400 520]);
    for s = 1:numel(iShow)
        i = iShow(s); subplot(1,numel(iShow),s); hold on;
        F = reshape(PhiGrid*snapBetaFull{io,i},size(GX));
        fvS = isosurface(GX,GY,GZ,F,1);
        fvS.vertices = fvS.vertices+cFull;
        pa = patch(fvS);
        set(pa,'FaceColor',[0.1 0.65 0.35],'EdgeColor','none','FaceAlpha',0.45);
        plot3(Pvis(1:4:end,1),Pvis(1:4:end,2),Pvis(1:4:end,3),'.', ...
            'Color',[0 0.3 0.8],'MarkerSize',3);
        plot3(Phid(1:4:end,1),Phid(1:4:end,2),Phid(1:4:end,3),'.', ...
            'Color',[0.6 0.6 0.6],'MarkerSize',3);
        plot3(cFull(1),cFull(2),cFull(3),'mx','MarkerSize',11,'LineWidth',2);
        axis equal; grid on; view(135,20); camlight; lighting gouraud;
        xlim([-1.6 1.6]); ylim([-1.6 1.6]); zlim([-1.6 1.6]);
        if lambdas(i)==0, tt='OLS (\lambda=0)'; else, tt=sprintf('\\lambda=%.0e',lambdas(i)); end
        title(sprintf('%s  CD_{occ}=%.3f',tt,CDOF(io,1,i)));
    end
    sgtitle(sprintf('occ=%d%%, full visible N=%d, view 1 — green: fit, magenta: center', ...
        round(100*occ_list(io)),size(Pvis,1)),'FontSize',12);
    saveFigBoth(fig,figDir,sprintf('fig3_full_o%02d',round(100*occ_list(io))));
    close(fig);
end

%% ---- Figure 4: ridge 이득 vs occlusion (clean-GT occ/all) ----
% 이득 = view별 (CD@λ=0 − min_λ CD), 표시 k는 kFocus
benOccV = nan(nOcc, nViews); benAllV = nan(nOcc, nViews);
for io = 1:nOcc
    for iv = 1:nViews
        cOcc = squeeze(CDO(io,iv,kFocusIdx,:));
        cAll = squeeze(CDA(io,iv,kFocusIdx,:));
        benOccV(io,iv) = cOcc(1) - min(cOcc);
        benAllV(io,iv) = cAll(1) - min(cAll);
    end
end
fig = figure('Visible','off','Position',[60 60 720 520]); hold on; grid on;
mbO = mean(benOccV,2,'omitnan'); sbO = std(benOccV,0,2,'omitnan');
mbA = mean(benAllV,2,'omitnan'); sbA = std(benAllV,0,2,'omitnan');
errorbar(100*occ_list(:),mbO(:),sbO(:),'o-','LineWidth',1.8, ...
    'Color',[0 0.3 0.8],'MarkerFaceColor',[0 0.3 0.8],'CapSize',8,'DisplayName','hidden region');
errorbar(100*occ_list(:),mbA(:),sbA(:),'s--','LineWidth',1.8, ...
    'Color',[0.8 0.25 0.1],'MarkerFaceColor',[0.8 0.25 0.1],'CapSize',8,'DisplayName','all surface');
xlabel('occlusion ratio [%]'); ylabel('ridge benefit: CD(OLS) - min_\lambda CD');
title(sprintf('%s ridge benefit, k=%.1f, mean\\pmstd / %d views',shapeName,kFocus,nViews));
xlim([max(0,100*min(occ_list)-5), min(100,100*max(occ_list)+5)]); legend('Location','best');
saveFigBoth(fig, figDir, 'fig4_benefit');

%% ---- Figure 5: df vs λ (focus k의 뷰 표준편차 + 전체 visible baseline) ----
fig = figure('Visible','off','Position',[50 50 1400 460]);
for p = 1:nPanel
    io = panelIdx(p);
    subplot(1, nPanel, p); hold on; grid on;
    dfM = squeeze(mDF(io,kFocusIdx,:)); dfS = squeeze(sDF(io,kFocusIdx,:));
    dfM = dfM(:); dfS = dfS(:);
    fill([lamPlot fliplr(lamPlot)], ...
        [max(dfM-dfS,0).' fliplr(min(dfM+dfS,30).')], [0.1 0.4 0.85], ...
        'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    if kFocus == 1.0, useI = 2:nL; else, useI = 1:nL; end
    penM = squeeze(mPEN(io,kFocusIdx,useI)); mseM = squeeze(mMSE(io,kFocusIdx,useI));
    iElb = useI(mengerElbow(log(penM.'),log(mseM.')));
    plot(lamPlot,dfM,'o-','Color',[0.1 0.4 0.85],'LineWidth',1.8, ...
        'DisplayName',sprintf('subset k=%.1f (elbow %.0e)',kFocus,lambdas(iElb)));
    plot(lamPlot(iElb),dfM(iElb),'s','Color',[0.1 0.4 0.85], ...
        'MarkerFaceColor',[0.1 0.4 0.85],'MarkerSize',9,'HandleVisibility','off');

    dfF = mDFF(io,:).'; dfFS = sDFF(io,:).';
    fill([lamPlot fliplr(lamPlot)], ...
        [max(dfF-dfFS,0).' fliplr(min(dfF+dfFS,30).')], [0.35 0.35 0.35], ...
        'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');
    iElbF = mengerElbow(log(mPENF(io,:)),log(mMSEF(io,:)));
    plot(lamPlot,dfF,'k--','LineWidth',2.0, ...
        'DisplayName',sprintf('all visible (elbow %.0e)',lambdas(iElbF)));
    plot(lamPlot(iElbF),dfF(iElbF),'ks','MarkerFaceColor','k', ...
        'MarkerSize',9,'HandleVisibility','off');

    set(gca, 'XScale', 'log');
    ylim([0 30]);
    xlabel('\lambda (OLS at 3e-5)'); ylabel('df(\lambda) = tr(H_\lambda)');
    title(sprintf('occ=%d%% (vis %d%%)', ...
        round(100*occ_list(io)), round(100*(1-occ_list(io)))));
    legend('Location','best','FontSize',7);
end
sgtitle(sprintf('%s effective dof — band: view std, square: L-curve elbow',shapeName));
saveFigBoth(fig, figDir, 'fig5_df');

%% ---- Figure 5b: 모든 subset k의 df 곡선 (표준편차 없는 보조 비교) ----
fig = figure('Visible','off','Position',[50 50 1400 460]);
for p = 1:nPanel
    io = panelIdx(p); subplot(1,nPanel,p); hold on; grid on;
    for kk = 1:nK
        plot(lamPlot,squeeze(mDF(io,kk,:)),'o-','Color',colsK(kk,:), ...
            'LineWidth',1.5,'DisplayName',sprintf('k=%.1f',k_list(kk)));
    end
    plot(lamPlot,mDFF(io,:),'k--','LineWidth',2,'DisplayName','all visible');
    set(gca,'XScale','log'); ylim([0 30]);
    xlabel('\lambda (OLS at 3e-5)'); ylabel('df(\lambda)');
    title(sprintf('occ=%d%%',round(100*occ_list(io)))); legend('Location','best','FontSize',7);
end
sgtitle(sprintf('%s effective dof — all subset sizes',shapeName));
saveFigBoth(fig,figDir,'fig5b_df_all_k');

%% ---- 콘솔 요약 ----
fprintf('\n=== %s summary (Sobolev, focus k=%.1f, mean over %d views) ===\n', ...
    shapeName, kFocus, nViews);
for io = 1:nOcc
    cO = squeeze(mCDO(io,kFocusIdx,:)); cA = squeeze(mCDA(io,kFocusIdx,:));
    if kFocus == 1.0, useI = 2:nL; else, useI = 1:nL; end
    penM = squeeze(mPEN(io,kFocusIdx,useI)); mseM = squeeze(mMSE(io,kFocusIdx,useI));
    iElb = useI(mengerElbow(log(penM.'),log(mseM.')));
    iElbF = mengerElbow(log(mPENF(io,:)),log(mMSEF(io,:)));
    [~,iOcc] = min(cO); [~,iAll] = min(cA);
    dfSub = mDF(io,kFocusIdx,iElb); dfFull = mDFF(io,iElbF);
    fprintf(['occ=%2d%% | CD_occ %.3f -> %.3f @λ=%.0e | CD_all %.3f -> %.3f @λ=%.0e | ' ...
             'df*_sub %.2f, df*_full %.2f, gap %+.2f\n'], ...
        round(100*occ_list(io)),cO(1),cO(iOcc),lambdas(iOcc), ...
        cA(1),cA(iAll),lambdas(iAll),dfSub,dfFull,dfSub-dfFull);
end
fprintf('empty hidden fitted surfaces: subset %d/%d, all-visible %d/%d (bbox-diagonal penalty)\n', ...
    nnz(OSF==0),numel(OSF),nnz(OSFF==0),numel(OSFF));
fprintf('saved figures + CSV to %s\n', figDir);

%% ---- local functions ----
% makeShapePoints / meshSurfacePoints / symmetricMetrics 는 utils/ 로 승격됨.
% Sec 3.2 가 같은 함수를 호출해야 cd_occ/cd_all 이 동일 정의임이 보장되고,
% 3.1→3.2 손실분해 사다리 비교가 성립한다. 여기에 다시 로컬로 두지 말 것
% (MATLAB 로컬 함수가 경로 함수를 가려 조용히 갈라진다).

function saveFigBoth(fig, figDir, name)
    exportgraphics(fig, fullfile(figDir, [name '.png']), 'Resolution', 150);
    savefig(fig, fullfile(figDir, [name '.fig']));
end

function iElb = mengerElbow(x, y)
    x = (x - min(x)) / max(max(x) - min(x), eps);
    y = (y - min(y)) / max(max(y) - min(y), eps);
    keep = 1;
    for i = 2:numel(x)
        if hypot(x(i)-x(keep(end)), y(i)-y(keep(end))) > 0.02
            keep(end+1) = i; %#ok<AGROW>
        end
    end
    n = numel(keep); kap = zeros(n,1);
    for j = 2:n-1
        a=[x(keep(j-1)) y(keep(j-1))]; b=[x(keep(j)) y(keep(j))]; c=[x(keep(j+1)) y(keep(j+1))];
        ar2 = abs((b(1)-a(1))*(c(2)-a(2)) - (c(1)-a(1))*(b(2)-a(2)));
        d = norm(b-a)*norm(c-b)*norm(c-a);
        if d>0, kap(j) = 2*ar2/d; end
    end
    [~, jMax] = max(kap);
    iElb = keep(jMax);
end
