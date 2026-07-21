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

% ---- 출력 폴더: 기존 규약대로 날짜 스탬프 ----
outDir   = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'test_ridge');
if ~isfolder(outDir), mkdir(outDir); end
runStamp = datestr(now, 'yyyymmdd_HHMM');
figDir   = fullfile(outDir, runStamp);
mkdir(figDir);
fprintf('그림 폴더 → %s\n', figDir);

%% ---- 기저 ----
order  = 6;
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
occ_list  = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70];   % 가려지는(hidden) 비율 (촘촘)
nOcc      = numel(occ_list);
snapOcc   = [0.10, 0.30, 0.50, 0.70]; % fig3 등위면 스냅샷을 낼 occ (10% 포함)
panelOcc  = [0.30, 0.50, 0.70];      % fig1/2/5 subplot으로 그릴 occ (fig4만 전체 7단계 사용)
nViews    = 10;
k_list    = [1.0, 1.2, 1.4];
nK        = numel(k_list);
lambdas   = [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1];
nL        = numel(lambdas);
nDrawScal = 10;                      % 스칼라 지표용 서브샘플 반복
hdPctl    = 95;                      % 로버스트 Hausdorff 백분위수

% 베이스 점군 1회 생성. 같은 10개 회전/노이즈 점군을 모든 occlusion에서 공유한다.
Pbase = generateRandomPointsOnHexagonPrism(totalN);
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
HDO = nan(nOcc, nViews, nK, nL);  DFm = nan(nOcc, nViews, nK, nL);
% fig3 스냅샷용: occ=0.50/0.70, view 1, 모든 k와 lambda의 β
snapBeta = cell(nOcc, nK, nL);
snapVis  = cell(nOcc, 1);
snapHid  = cell(nOcc, 1);
snapCen  = cell(nOcc, nK);   % 스냅 드로우 무게중심(등위면 원프레임 복원용)

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
        Phid   = Pnoisy(~vis, :);

        NvisAll = size(Pvis, 1);

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

                % --- 거리 지표: rep1 드로우의 등위면 (센터링 프레임 → 원프레임 복원) ---
                Phi_s = PhiPool{1};
                bCD = (Phi_s.'*Phi_s + D) \ (Phi_s.'*ones(Nfit,1));
                F   = reshape(PhiGrid*bCD, size(GX));
                fv  = isosurface(GX, GY, GZ, F, 1);
                V   = fv.vertices + cCD;   % 등위면 정점을 원프레임으로 복원(+cCD)
                if size(V,1) > 15000, V = V(randperm(size(V,1),15000), :); end
                if ~isempty(V)
                    % CD와 HD를 방향별 nnDist 1회씩으로 동시 도출
                    d1 = nnDist(Phid, V);      % hidden → surface
                    d2 = nnDist(V, Phid);      % surface → hidden
                    CDO(io,iv,kk,i) = 0.5*(mean(d1) + mean(d2));
                    HDO(io,iv,kk,i) = max(prctile(d1,hdPctl), prctile(d2,hdPctl));
                end

                Trows(end+1,:) = {occ, iv, k_list(kk), lam, ...
                    MSE(io,iv,kk,i), PEN(io,iv,kk,i), GCVv(io,iv,kk,i), ...
                    DFm(io,iv,kk,i), CDO(io,iv,kk,i), HDO(io,iv,kk,i)}; %#ok<AGROW>

                if any(abs(occ - snapOcc) < 1e-9) && iv == 1
                    snapBeta{io,kk,i} = bCD;
                    snapCen{io,kk}    = cCD;
                    if isempty(snapVis{io})
                        snapVis{io} = Pvis;
                        snapHid{io} = Phid;
                    end
                end
            end
        end
        fprintf('occ=%.2f view=%2d/%d  (vis %d / hid %d)\n', ...
            occ, iv, nViews, size(Pvis,1), size(Phid,1));
    end
end

% 뷰 집계
mCDO = squeeze(mean(CDO, 2, 'omitnan'));  sCDO = squeeze(std(CDO, 0, 2, 'omitnan'));
mHDO = squeeze(mean(HDO, 2, 'omitnan'));  sHDO = squeeze(std(HDO, 0, 2, 'omitnan'));
mMSE = squeeze(mean(MSE, 2, 'omitnan'));  mPEN = squeeze(mean(PEN, 2, 'omitnan'));
mGCV = squeeze(mean(GCVv,2, 'omitnan'));  % [occ,k,λ]
mDF  = squeeze(mean(DFm, 2, 'omitnan'));  % [occ,k,λ]

%% ---- CSV ----
T = cell2table(Trows, 'VariableNames', ...
    {'occ','view','k','lambda','mse','pen','gcv','df','cd_occ','hd_occ'});
writetable(T, fullfile(figDir, 'lcurve_demo_results.csv'));

%% ---- 그림 축 준비 ----
lamPlot = lambdas; lamPlot(1) = 3e-5;   % λ=0 로그축 표시용
colsK   = lines(nK);
kMid    = find(k_list == 1.2);
panelIdx = find(ismember(round(occ_list*100), round(panelOcc*100)));  % fig1/2/5용 occ 인덱스
nPanel   = numel(panelIdx);

%% ---- Figure 1: L-curve (occ별 subplot, 모든 k, 뷰 평균+스프레드) ----
% k=1.0(N=nT=28)은 λ=0에서 보간→MSE≈0 이상치라 λ=0 점을 제외하고 그린다.
% y축 하한은 1e-2로 고정(의미 있는 점 min≈2.8e-2은 모두 유지, 보간 이상치만 배제).
yTop = max(mMSE(:)) * 1.5;
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
        if kk == kMid
            for i = [2 4 6 8]
                text(mPEN(io,kk,i)*1.06, mMSE(io,kk,i), sprintf('%.0e',lambdas(i)), ...
                    'Color', colsK(kk,:), 'FontSize', 7);
            end
        end
    end
    set(gca,'XScale','log','YScale','log');
    ylim([1e-2, yTop]);
    xlabel('\beta^T D \beta'); ylabel('algebraic MSE');
    title(sprintf('occ=%d%% (vis %d%%)', ...
        round(100*occ_list(io)), round(100*(1-occ_list(io)))));
    legend('Location','best','FontSize',7);
end
sgtitle('L-curve (Sobolev) — color: k, thin: 10 paired views, bold: mean, square: elbow  (k=1.0: λ=0 보간점 제외)');
saveFigBoth(fig, figDir, 'fig1_lcurve');

%% ---- Figure 2: λ vs CD_occ / HD_occ (2×nPanel, mean±std, k별) ----
fig = figure('Visible','off','Position',[40 40 1400 760]);
metNames = {'CD_{occ}','HD_{occ} (95pct)'};
mMet = {mCDO, mHDO};  sMet = {sCDO, sHDO};
for mrow = 1:2
    for p = 1:nPanel
        io = panelIdx(p);
        subplot(2, nPanel, (mrow-1)*nPanel + p); hold on; grid on;
        for kk = 1:nK
            mu = squeeze(mMet{mrow}(io,kk,:)); sd = squeeze(sMet{mrow}(io,kk,:));
            fill([lamPlot fliplr(lamPlot)], [(mu-sd).' fliplr((mu+sd).')], ...
                colsK(kk,:), 'FaceAlpha', 0.12, 'EdgeColor','none', 'HandleVisibility','off');
            plot(lamPlot, mu, 'o-', 'Color', colsK(kk,:), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('k=%.1f', k_list(kk)));
            [~,imin] = min(mu);
            plot(lamPlot(imin), mu(imin), 'p', 'Color', colsK(kk,:), ...
                'MarkerFaceColor', colsK(kk,:), 'MarkerSize', 10, 'HandleVisibility','off');
        end
        set(gca,'XScale','log','YScale','log');
        if mrow==2, xlabel('\lambda (OLS at 3e-5)'); end
        ylabel(metNames{mrow});
        if mrow==1, title(sprintf('occ=%d%%', round(100*occ_list(io)))); end
        if p==1 && mrow==1, legend('Location','northwest','FontSize',8); end
    end
end
sgtitle('occluded-region error vs \lambda — mean\pmstd over 10 views (star = min)');
saveFigBoth(fig, figDir, 'fig2_cd_hd');

%% ---- Figure 3: 등위면 스냅샷 (snapOcc=30/50/70%, 모든 k, view1) ----
% 각 그림은 OLS, lambda=1e-3, lambda=1e-1의 3-panel 비교다.
iShow = [find(lambdas==0), find(lambdas==1e-3), find(lambdas==1e-1)];
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
            plot3(Pvis(1:4:end,1),Pvis(1:4:end,2),Pvis(1:4:end,3),'.', ...
                'Color',[0 0.3 0.8],'MarkerSize',3);
            plot3(Phid(1:4:end,1),Phid(1:4:end,2),Phid(1:4:end,3),'.', ...
                'Color',[0.6 0.6 0.6],'MarkerSize',3);
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
            'blue: visible, gray: hidden'], ...
            round(100*occ_list(io)), k_list(kk), Nfit));
        figName = sprintf('fig3_o%02d_k%02d', ...
            round(100*occ_list(io)), round(10*k_list(kk)));
        saveFigBoth(fig, figDir, figName);
        close(fig);
    end
end

%% ---- Figure 4: ridge 이득 vs occlusion (10뷰 에러바) ----
% 이득 = view별 (CD_occ@λ=0 − min_λ CD_occ), k=1.2
benV = nan(nOcc, nViews);
for io = 1:nOcc
    for iv = 1:nViews
        c = squeeze(CDO(io,iv,kMid,:));
        benV(io,iv) = c(1) - min(c);
    end
end
fig = figure('Visible','off','Position',[60 60 720 520]); hold on; grid on;
mb = mean(benV,2,'omitnan'); sb = std(benV,0,2,'omitnan');
errorbar(100*occ_list(:), mb(:), sb(:), 'o-', 'LineWidth', 1.8, ...
    'Color', [0 0.3 0.8], 'MarkerFaceColor', [0 0.3 0.8], 'CapSize', 10);
xlabel('occlusion ratio [%]'); ylabel('ridge benefit:  CD_{occ}(OLS) − min_\lambda CD_{occ}');
title('Ridge benefit vs occlusion (Sobolev, k=1.2, mean\pmstd / 10 views)');
xlim([5 75]);
saveFigBoth(fig, figDir, 'fig4_benefit');

%% ---- Figure 5: df vs λ (유효 자유도 진단) ----
% df(λ)=tr(H_λ). occ별 subplot, k별 라인. L-curve elbow λ에 마커.
% df는 적합/L-curve를 바꾸지 않는 사후 진단값(변경 없이 새 축만 추가).
fig = figure('Visible','off','Position',[50 50 1400 460]);
for p = 1:nPanel
    io = panelIdx(p);
    subplot(1, nPanel, p); hold on; grid on;
    for kk = 1:nK
        dfM  = squeeze(mDF(io,kk,:));
        % elbow는 fig1과 동일 규약: k=1.0은 λ=0(보간 이상치) 제외
        if k_list(kk) == 1.0, useI = 2:nL; else, useI = 1:nL; end
        penM = squeeze(mPEN(io,kk,useI)); mseM = squeeze(mMSE(io,kk,useI));
        iElb = useI(mengerElbow(log(penM.'), log(mseM.')));   % 전체 인덱스로 복원
        plot(lamPlot, dfM, 'o-', 'Color', colsK(kk,:), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('k=%.1f (elbow %.0e)', k_list(kk), lambdas(iElb)));
        plot(lamPlot(iElb), dfM(iElb), 's', 'Color', colsK(kk,:), ...
            'MarkerFaceColor', colsK(kk,:), 'MarkerSize', 9, 'HandleVisibility','off');
    end
    set(gca, 'XScale', 'log');
    ylim([0 30]);                         % 모든 occ 패널에서 df 축 범위 통일
    xlabel('\lambda (OLS at 3e-5)'); ylabel('df(\lambda) = tr(H_\lambda)');
    title(sprintf('occ=%d%% (vis %d%%)', ...
        round(100*occ_list(io)), round(100*(1-occ_list(io)))));
    legend('Location','best','FontSize',7);
end
sgtitle('effective dof  df(\lambda)  vs  \lambda — color: k, square: L-curve elbow');
saveFigBoth(fig, figDir, 'fig5_df');

%% ---- 콘솔 요약 ----
fprintf('\n=== summary (Sobolev, k=1.2, mean over %d views) ===\n', nViews);
for io = 1:nOcc
    cM = squeeze(mCDO(io,kMid,:)); hM = squeeze(mHDO(io,kMid,:));
    penM = squeeze(mPEN(io,kMid,:)); mseM = squeeze(mMSE(io,kMid,:));
    iElb = mengerElbow(log(penM.'), log(mseM.'));
    [~,iCD] = min(cM); [~,iHD] = min(hM);
    fprintf(['occ=%2d%% | CD OLS %.3f -> min %.3f @λ=%.0e | HD min %.3f @λ=%.0e | ' ...
             'elbow λ=%.0e (CD %.3f) | benefit %.3f\n'], ...
        round(100*occ_list(io)), cM(1), cM(iCD), lambdas(iCD), hM(iHD), lambdas(iHD), ...
        lambdas(iElb), cM(iElb), cM(1)-cM(iCD));
end
fprintf('saved figures + CSV to %s\n', figDir);

%% ---- local functions ----
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
