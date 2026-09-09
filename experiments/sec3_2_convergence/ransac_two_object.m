%% ransac_two_object.m — 2물체 구조적 아웃라이어 RANSAC (Sec 3.2)
%
% 타깃 A(다수) + 방해물 B(소수, 구조적 아웃라이어)를 한 장면에 배치.
% RANSAC은 A에 lock-on 해야 하며 B는 인라이어처럼 보이는 어려운 오염이다.
% gtMask=A. ε=n_B/(n_A+n_B). CD/HD GT = A의 가림 없는 clean 전체 표면.
%
% 지표는 ransac_outlier_sweep.m 과 동일 체계(같은 utils 함수) + 2물체 고유:
%   lock_on      — 인라이어의 과반이 A인가. '경계 혼동' vs '방해물에 물림' 구분용.
%   e_par/e_off  — 중심 편향을 가림 방향(월드 z)과 방해물 방향(뷰프레임 x)으로 분해.
%   fp_lev       — 오분류 점의 GT 표면까지 거리.
%
% 2물체 occlusion 설계 (switch):
%   'perobj' (주) : A,B 각각 z-상위 (1-occ) 남기고 합침. 타깃 가시율이 배치와 무관하게 occ.
%   'global'(보조): 합친 뒤 전역 z-상위 (1-occ). 먼 물체가 더 가려진다(상호가림).
% 두 도형은 동일 view 회전(단일 시점) 적용. 오프셋 d로 A,B 무게중심 분리.
%
% 실행: MATLAB R2025a.  override: quickMode=true; seedOffset=100; run(...)

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
if ~exist('quickMode','var')  || isempty(quickMode),  quickMode  = false; end
if ~exist('seedOffset','var') || isempty(seedOffset), seedOffset = 0;     end

%% ---- 공통 파라미터 ----
totalN   = 4000; noiseSig = 0.007; occ = 0.30;
conf = 0.95; kMult = 1.4; tau = 0.45;
failFrac = 0.5; gapMin = 1.15;
offMap = struct('overlap',0.6, 'adjacent',1.6, 'separate',2.8);  % 무게중심 오프셋 거리

if quickMode
    nSeed=3; eps_axis=[0.2 0.3]; off_axis={'adjacent'};
    scheme_axis={'perobj'}; pair_axis={{'hexagon','cube'}}; maxIterV=300;
else
    nSeed=20; eps_axis=[0.1 0.2 0.3 0.4]; off_axis={'overlap','adjacent','separate'};
    scheme_axis={'perobj','global'};
    pair_axis={{'hexagon','cube'},{'hexagon','cylinder'},{'cube','cylinder'}};
    maxIterV=1000;
end
baseEps=0.3; baseOff='adjacent'; baseScheme='perobj'; basePair={'hexagon','cube'};
baseLam=1e-2; lam6=3e-1;   % 4차 / 6차 기준 λ

%% ---- 기저 ----
order4=4; TermsC4=homogeneFischerTerms(order4); [Funcs4,Grads4,nT4]=makeFuncsGradsStack(TermsC4);
order6=6; TermsC6=homogeneFischerTerms(order6); [Funcs6,~,nT6]=makeFuncsGradsStack(TermsC6);
Nfit=round(kMult*nT4);
gLin=linspace(-2.2,2.2,55); [GX,GY,GZ]=meshgrid(gLin,gLin,gLin);
PhiGrid4=calculateFourthOrder([GX(:),GY(:),GZ(:)],Funcs4);
PhiGrid6=calculateFourthOrder([GX(:),GY(:),GZ(:)],Funcs6);
nSurf=12000; hdPct=95;

ransacPar=struct('maxIter',maxIterV,'conf',conf,'thresh',tau, ...
    'minInlierRatio',0.65,'updateThresh',0,'locIters',4, ...
    'momentum',0.66,'damping',0.85,'reg',1e-6,'k',kMult,'lambda',baseLam);
ransacPar.mc=struct('on',true,'saveVarName','','saveMatFile','','metric','basic');
ransacPar.sim=struct('freezeW',false,'freezeIter',true,'freezeOmega',false);

%% ---- OFAT 셀 ----
cellList={};  % {pair, eps, off, scheme}
for e=eps_axis,   cellList(end+1,:)={basePair,e,baseOff,baseScheme}; end %#ok<*SAGROW>
for oi=1:numel(off_axis),    cellList(end+1,:)={basePair,baseEps,off_axis{oi},baseScheme}; end
for si=1:numel(scheme_axis), cellList(end+1,:)={basePair,baseEps,baseOff,scheme_axis{si}}; end
for pi2=1:numel(pair_axis),  cellList(end+1,:)={pair_axis{pi2},baseEps,baseOff,baseScheme}; end
key=cellfun(@(p,e,o,s) sprintf('%s+%s|%.2f|%s|%s',p{1},p{2},e,o,s), ...
    cellList(:,1),cellList(:,2),cellList(:,3),cellList(:,4),'uni',0);
[~,ia]=unique(key,'stable'); cellList=cellList(ia,:); nCell=size(cellList,1);

%% ---- 출력 ----
outDir=fullfile(fileparts(mfilename('fullpath')),'..','..','test_ridge');
stamp=sprintf('sec32_twoobj_%s%s%s',datestr(now,'yyyymmdd_HHMM'), ...
    tern(seedOffset>0,sprintf('_fold%d',seedOffset),''), tern(quickMode,'_quick',''));
figDir=fullfile(outDir,stamp); mkdir(figDir);
fprintf('출력 폴더 → %s\n',figDir);
fprintf('셀 %d개 × 시드 %d회 = 실행 %d회 (seedOffset=%d)\n', nCell,nSeed,nCell*nSeed,seedOffset);
% -batch 실행은 stdout 이 버퍼링되어 밖에서 진행을 못 본다. 셀마다 한 줄 append.
%   확인:  Get-Content test_ridge\<런폴더>\progress.log -Tail 5 -Wait
logPath=fullfile(figDir,'progress.log');
logLine(logPath,'START %s | cells=%d seeds=%d seedOffset=%d', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS'),nCell,nSeed,seedOffset);
tRun=tic;

%% ---- 실행 ----
Rows={};
for ci=1:nCell
    pair=cellList{ci,1}; epsv=cellList{ci,2}; offn=cellList{ci,3}; sch=cellList{ci,4};
    d=offMap.(offn); ransacPar.lambda=baseLam; ransacPar.maxIter=maxIterV;
    isBase=isequal(pair,basePair)&&abs(epsv-baseEps)<1e-9&&strcmp(offn,baseOff)&&strcmp(sch,baseScheme);
    for sd=1:nSeed
        rng(sd + seedOffset);
        % 타깃 A 점수(다수), 방해물 B 점수(소수)를 ε로 배분
        nA=round(totalN*(1-epsv)); nB=round(totalN*epsv);
        Rv=randomRotation();               % 단일 시점: 두 도형 공통 회전
        A0=makeShapePoints(pair{1},nA)*Rv.';
        B0=makeShapePoints(pair{2},nB)*Rv.';
        offDir=[1 0 0]*Rv.';               % 오프셋 방향(뷰 프레임 x)
        B0=B0 + d*offDir;                  % 방해물 이동
        Aclean=A0;                          % GT = A의 clean 전체(가림 전)
        bboxDiag=norm(max(Aclean,[],1)-min(Aclean,[],1));

        % --- occlusion (visA 마스크를 유지해 GT 3분할에 사용) ---
        switch sch
            case 'perobj'   % 각자 z-상위 (1-occ)
                zThrA = prctile(A0(:,3),100*occ);
                visA  = A0(:,3) >= zThrA;
                visB  = B0(:,3) >= prctile(B0(:,3),100*occ);
            case 'global'   % 합친 뒤 전역 z-상위 (1-occ)
                zThrA = prctile([A0(:,3);B0(:,3)],100*occ);
                visA  = A0(:,3) >= zThrA;
                visB  = B0(:,3) >= zThrA;
        end
        Avis=A0(visA,:); Bvis=B0(visB,:);
        Pin=[Avis;Bvis]+noiseSig*randn(size(Avis,1)+size(Bvis,1),3);
        gtMask=[true(size(Avis,1),1); false(size(Bvis,1),1)];
        ransacPar.mc.gtMask=gtMask;

        [~]=evalc(['[ResultDisp,Beta,inMask,Log] = ' ...
            'PoliNavigationSolver3_3_FischerRansac_MC(0,Pin,order4,nT4,Funcs4,Grads4,ransacPar);']);
        inMask=logical(inMask);
        TP=sum(gtMask&inMask); FP=sum(~gtMask&inMask); FN=sum(gtMask&~inMask);
        prec=TP/max(1,TP+FP); rec=TP/max(1,TP+FN); F1=2*prec*rec/max(1e-12,prec+rec);

        % --- 2물체 고유: 어느 물체에 물었는가 ---
        inlier_A_frac = TP / max(1,sum(inMask));
        lock_on = double(inlier_A_frac > 0.5);

        ps=[Log.postScore]'; ps(isnan(ps))=-inf; rb=cummax(ps);
        nConv=find(rb>=0.99*rb(end),1,'first'); if isempty(nConv), nConv=maxIterV; end

        % --- τ 게이트 진단 ---
        rGt = abs(calculateFourthOrder(Pin(gtMask,:) - ResultDisp, Funcs4)*Beta - 1);
        r_med_gt = median(rGt); gate_margin = r_med_gt/tau;

        % --- L3: 실제 RANSAC 기하 (A의 clean 전체 대비, 3분할) ---
        S3 = splitSurfaceMetrics(Beta, PhiGrid4,GX,GY,GZ, ResultDisp, nSurf, ...
                                 Aclean, visA, zThrA, hdPct);
        % --- L2 (oracle): 같은 솔버에 A의 가시점만 넣은 상한 ---
        % 1회 적합으로 대체하면 중심이 가시점 평균에 고정돼 상한이 오히려 나빠진다.
        Porc=Pin(gtMask,:);
        parOracle=ransacPar; parOracle.mc.gtMask=true(size(Porc,1),1);
        [~]=evalc(['[dispO,betaO,~,~] = PoliNavigationSolver3_3_FischerRansac_MC' ...
            '(0,Porc,order4,nT4,Funcs4,Grads4,parOracle);']);
        S2 = splitSurfaceMetrics(betaO, PhiGrid4,GX,GY,GZ, dispO, nSurf, ...
                                 Aclean, visA, zThrA, hdPct);
        % --- L2d: RANSAC 없이 1회만 적합. L2 - L2d = 중심 최적화 몫 ---
        cO=mean(Porc,1);
        bD=regressionFourthOrder(Porc-cO,Funcs4,baseLam,order4,1);
        Sd = splitSurfaceMetrics(bD, PhiGrid4,GX,GY,GZ, cO, nSurf, ...
                                 Aclean, visA, zThrA, hdPct);

        % --- 오분류의 기하 ---
        fpIdx=~gtMask&inMask;
        if any(fpIdx)
            dFP=nnDist(Pin(fpIdx,:),Aclean);
            fp_lev=median(dFP); fp_lev_p95=prctile(dFP,95);
        else
            fp_lev=0; fp_lev_p95=0;
        end
        fn_clump = spreadNN(Pin(gtMask&~inMask,:)) / spreadNN(Pin(gtMask,:));
        % 단일물체판과 달리 중심오차 기준을 안 쓴다. 2물체에서는 제대로 물고도 중심이
        % 방해물 쪽으로 끌리므로(e_off 최대 27%) 오검출이 난다. 검출기는 lock_on.
        failFlag = S3.empty || S3.cd_all > failFrac*bboxDiag;

        % --- 중심 편향: 가림 방향 / 방해물 방향 / 잔차 로 직교 분해 ---
        % offDir 은 월드 z 와 직교하지 않으므로 Gram-Schmidt 로 정규직교화한다.
        eVec = ResultDisp(:).' - mean(Aclean,1);
        u1 = [0 0 1];
        u2 = offDir - dot(offDir,u1)*u1;
        if norm(u2) > 1e-9, u2 = u2/norm(u2); else, u2 = [1 0 0]; end
        e_par = dot(eVec,u1)/bboxDiag;                     % 가림 방향
        e_off = dot(eVec,u2)/bboxDiag;                     % 방해물 방향(가림-직교 성분)
        e_res = norm(eVec - dot(eVec,u1)*u1 - dot(eVec,u2)*u2)/bboxDiag;

        % --- 장축 자세 (GT/적합 양쪽 게이트) ---
        [axGt,gapGt]=principalAxis(Aclean);
        gapFit=NaN; theta_axis=NaN; theta_axis_pca=NaN;
        if ~S3.empty
            [axFit,gapFit]=principalAxis(S3.V);
            if gapGt>gapMin
                theta_axis_pca=acosd(min(1,abs(dot(principalAxis(Avis),axGt))));
                if gapFit>gapMin
                    theta_axis=acosd(min(1,abs(dot(axFit,axGt))));
                end
            end
        end

        % --- 6차 post-hoc refit (기준셀만) ---
        cd_occ6=NaN; cd_all6=NaN;
        if isBase
            Pinl=Pin(inMask,:); cc=mean(Pinl,1);
            b6=regressionFourthOrder(Pinl-cc,Funcs6,lam6,order6,1);
            S6=splitSurfaceMetrics(b6,PhiGrid6,GX,GY,GZ,cc,nSurf,Aclean,visA,zThrA,hdPct);
            cd_occ6=S6.cd_occ; cd_all6=S6.cd_all;
        end

        Rows(end+1,:)={sprintf('%s+%s',pair{1},pair{2}),epsv,offn,sch,d,sd+seedOffset, ...
            size(Avis,1),size(Bvis,1),prec,rec,F1,nConv, ...
            inlier_A_frac,lock_on, ...
            S3.cd_occ,S3.hd_occ,S3.cd_vis,S3.hd_vis,S3.cd_all,S3.hd_all,S3.osf, ...
            S2.cd_occ,S2.cd_all,Sd.cd_occ,Sd.cd_all, ...
            fp_lev,fp_lev_p95,fn_clump,double(failFlag), ...
            r_med_gt,gate_margin, ...
            e_par,e_off,e_res,theta_axis,theta_axis_pca,gapGt,gapFit, ...
            cd_occ6,cd_all6};
    end
    el=toc(tRun); eta=el/ci*(nCell-ci);
    msg=sprintf('[%2d/%2d] %s+%s eps=%.2f off=%-8s sch=%-6s',ci,nCell,pair{1},pair{2},epsv,offn,sch);
    fprintf('  %s done\n',msg);
    logLine(logPath,'%s done | 경과 %s, 남은 예상 %s',msg,hms(el),hms(eta));
end
logLine(logPath,'LOOP DONE %s | 총 %s. 이제 CSV/그림 저장 중...', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS'),hms(toc(tRun)));

T=cell2table(Rows,'VariableNames',{'pair','eps','offset','scheme','d','seed', ...
    'n_A','n_B','precision','recall','F1','n_conv', ...
    'inlier_A_frac','lock_on', ...
    'cd_occ','hd_occ','cd_vis','hd_vis','cd_all','hd_all','osf', ...
    'cd_occ_oracle','cd_all_oracle','cd_occ_direct','cd_all_direct', ...
    'fp_lev','fp_lev_p95','fn_clump','fail','r_med_gt','gate_margin', ...
    'e_par','e_off','e_res','theta_axis','theta_axis_pca','axis_gap_gt','axis_gap_fit', ...
    'cd_occ6','cd_all6'});
writetable(T,fullfile(figDir,'ransac_twoobj_results.csv'));
fprintf('CSV 저장 → %s\n',fullfile(figDir,'ransac_twoobj_results.csv'));

%% ---- 그림 ----
makeFiguresTwo(T,figDir,basePair,baseOff,baseScheme,baseEps);
fprintf('그림 저장 → %s\n',figDir);

%% ============================ local functions ============================
% makeShapePoints / randomRotation / meshSurfacePoints 는 utils/ 에 있다.
function s=tern(c,a,b), if c,s=a; else,s=b; end, end

function logLine(path,fmt,varargin)
% 한 줄 append 후 즉시 닫는다.
    fid=fopen(path,'a'); if fid<0, return; end
    fprintf(fid,[fmt '\n'],varargin{:}); fclose(fid);
end

function s=hms(sec)
    if ~isfinite(sec), s='--'; return; end
    h=floor(sec/3600); m=floor(mod(sec,3600)/60);
    if h>0, s=sprintf('%dh%02dm',h,m); else, s=sprintf('%dm%02ds',m,floor(mod(sec,60))); end
end

function makeFiguresTwo(T,figDir,basePair,baseOff,baseScheme,baseEps)
    bp=sprintf('%s+%s',basePair{1},basePair{2});
    med    = @(mask,var) median(T.(var)(mask & T.fail==0),'omitnan');
    rate   = @(mask,var) mean(T.(var)(mask));
    byOff  = @(e,o) strcmp(T.pair,bp)&abs(T.eps-e)<1e-9&strcmp(T.offset,o)&strcmp(T.scheme,baseScheme);
    epss=unique(T.eps)'; offs=unique(T.offset); schs=unique(T.scheme); prs=unique(T.pair);

    % fig1: F1 vs eps (offset별 / 스킴 / 도형쌍)
    fig=figure('Visible','off','Position',[50 50 1300 420]);
    subplot(1,3,1); hold on; grid on;
    for oi=1:numel(offs)
        y=arrayfun(@(e) med(byOff(e,offs{oi}),'F1'),epss);
        if all(isnan(y)),continue; end
        plot(100*epss,y,'-o','LineWidth',1.6,'DisplayName',offs{oi});
    end
    xlabel('\epsilon [%] (방해물 비율)'); ylabel('F1^{inlier}'); title('F1 vs \epsilon (offset별)');
    legend('Location','southwest','FontSize',7); ylim([0 1]);
    subplot(1,3,2); hold on; grid on;
    for si=1:numel(schs)
        m=@(e) strcmp(T.pair,bp)&abs(T.eps-e)<1e-9&strcmp(T.offset,baseOff)&strcmp(T.scheme,schs{si});
        y=arrayfun(@(e) med(m(e),'F1'),epss);
        if all(isnan(y)),continue; end
        plot(100*epss,y,'-o','LineWidth',1.6,'DisplayName',schs{si});
    end
    xlabel('\epsilon [%]'); ylabel('F1^{inlier}'); title('가림 스킴 비교'); legend('Location','southwest','FontSize',7); ylim([0 1]);
    subplot(1,3,3); hold on; grid on;
    for pj=1:numel(prs)
        m=@(e) strcmp(T.pair,prs{pj})&abs(T.eps-e)<1e-9&strcmp(T.offset,baseOff)&strcmp(T.scheme,baseScheme);
        y=arrayfun(@(e) med(m(e),'F1'),epss);
        if all(isnan(y)),continue; end
        plot(100*epss,y,'-o','LineWidth',1.6,'DisplayName',prs{pj});
    end
    xlabel('\epsilon [%]'); ylabel('F1^{inlier}'); title('도형쌍별'); legend('Location','southwest','FontSize',7); ylim([0 1]);
    saveFig(fig,figDir,'fig1_twoobj_F1');

    % fig2: CD 3분할 + 손실 사다리 (offset별)
    fig=figure('Visible','off','Position',[60 60 1300 420]);
    subplot(1,3,1); hold on; grid on;
    for oi=1:numel(offs)
        y=arrayfun(@(e) med(byOff(e,offs{oi}),'cd_all'),epss);
        if all(isnan(y)),continue; end
        plot(100*epss,y,'-o','LineWidth',1.6,'DisplayName',offs{oi});
    end
    xlabel('\epsilon [%]'); ylabel('CD\_all (타깃 A 대비)'); title('2물체 Chamfer vs \epsilon');
    legend('Location','northwest','FontSize',7);
    subplot(1,3,2); hold on; grid on;
    y3=arrayfun(@(e) med(byOff(e,baseOff),'cd_all'),epss);
    y2=arrayfun(@(e) med(byOff(e,baseOff),'cd_all_oracle'),epss);
    yd=arrayfun(@(e) med(byOff(e,baseOff),'cd_all_direct'),epss);
    fill([100*epss fliplr(100*epss)],[y2 fliplr(y3)],[0.85 0.4 0.4], ...
         'FaceAlpha',0.18,'EdgeColor','none','DisplayName','아웃라이어를 못 걸러서 생긴 손실');
    plot(100*epss,yd,'-^','LineWidth',1.4,'Color',[.5 .5 .5],'DisplayName','RANSAC 없이 1회 적합');
    plot(100*epss,y2,'-s','LineWidth',1.8,'DisplayName','정답만 넣었을 때 (도달 가능 상한)');
    plot(100*epss,y3,'-o','LineWidth',1.8,'DisplayName','실제 결과');
    xlabel('\epsilon [%]'); ylabel('CD\_all');
    title('손실 분해 (회색과의 차이 = 중심을 맞춰서 좋아진 양)'); legend('Location','northwest','FontSize',7);
    subplot(1,3,3); hold on; grid on;
    for oi=1:numel(offs)
        y=arrayfun(@(e) med(byOff(e,offs{oi}),'cd_occ'),epss);
        if all(isnan(y)),continue; end
        plot(100*epss,y,'-o','LineWidth',1.6,'DisplayName',offs{oi});
    end
    xlabel('\epsilon [%]'); ylabel('CD\_occ (가려진 영역)'); title('가려진 영역 복원'); legend('Location','northwest','FontSize',7);
    saveFig(fig,figDir,'fig2_twoobj_cd');

    % fig3: offset 축의 fp_lev / CD
    fig=figure('Visible','off','Position',[70 70 1300 420]);
    offOrder = {'overlap','adjacent','separate'};
    offOrder = offOrder(ismember(offOrder, offs));
    xo = 1:numel(offOrder);
    getv=@(o,v) med(strcmp(T.pair,bp)&abs(T.eps-baseEps)<1e-9&strcmp(T.offset,o)&strcmp(T.scheme,baseScheme),v);
    subplot(1,3,1); hold on; grid on;
    yyaxis left;  plot(xo,cellfun(@(o) getv(o,'F1'),offOrder),'-o','LineWidth',2); ylabel('F1^{inlier}'); ylim([0 1]);
    yyaxis right; plot(xo,cellfun(@(o) getv(o,'fp_lev'),offOrder),'-s','LineWidth',2); ylabel('fp\_lev');
    set(gca,'XTick',xo,'XTickLabel',offOrder); title('F1은 평평, fp\_lev는 상승');
    subplot(1,3,2); hold on; grid on;
    plot(xo,cellfun(@(o) getv(o,'cd_all'),offOrder),'-o','LineWidth',2,'DisplayName','CD\_all');
    plot(xo,cellfun(@(o) getv(o,'fp_lev_p95'),offOrder),'-s','LineWidth',2,'DisplayName','fp\_lev p95');
    set(gca,'XTick',xo,'XTickLabel',offOrder); title('잘못 넣은 점의 거리 → CD'); legend('Location','northwest','FontSize',7);
    subplot(1,3,3); hold on; grid on;
    lo=cellfun(@(o) rate(strcmp(T.pair,bp)&abs(T.eps-baseEps)<1e-9&strcmp(T.offset,o)&strcmp(T.scheme,baseScheme),'lock_on'),offOrder);
    fr=cellfun(@(o) rate(strcmp(T.pair,bp)&abs(T.eps-baseEps)<1e-9&strcmp(T.offset,o)&strcmp(T.scheme,baseScheme),'fail'),offOrder);
    bar(xo,[lo(:) fr(:)]); set(gca,'XTick',xo,'XTickLabel',offOrder);
    legend({'lock-on 성공률','failRate'},'Location','best','FontSize',7); ylim([0 1]);
    title('타깃 A 에 물었는가');
    saveFig(fig,figDir,'fig3_twoobj_mechanism');

    % fig4: 중심 편향의 방향 분해 (가림 vs 방해물)
    fig=figure('Visible','off','Position',[80 80 1000 420]);
    subplot(1,2,1); hold on; grid on;
    plot(xo,cellfun(@(o) getv(o,'e_par'),offOrder),'-o','LineWidth',2,'DisplayName','e_{par} (가림 방향)');
    plot(xo,cellfun(@(o) getv(o,'e_off'),offOrder),'-s','LineWidth',2,'DisplayName','e_{off} (방해물 방향)');
    plot(xo,cellfun(@(o) getv(o,'e_res'),offOrder),'-^','LineWidth',1.2,'DisplayName','e_{res} (잔차)');
    yline(0,'k--'); set(gca,'XTick',xo,'XTickLabel',offOrder);
    ylabel('normalized bias'); title('중심 편향 방향 분해'); legend('Location','best','FontSize',7);
    subplot(1,2,2); hold on; grid on;
    plot(100*epss,arrayfun(@(e) med(byOff(e,baseOff),'e_par'),epss),'-o','LineWidth',1.8,'DisplayName','e_{par}');
    plot(100*epss,arrayfun(@(e) med(byOff(e,baseOff),'e_off'),epss),'-s','LineWidth',1.8,'DisplayName','e_{off}');
    yline(0,'k--'); xlabel('\epsilon [%]'); ylabel('normalized bias');
    title('오염량에 따른 편향'); legend('Location','best','FontSize',7);
    saveFig(fig,figDir,'fig4_twoobj_bias');
end

function saveFig(fig,figDir,name)
    exportgraphics(fig,fullfile(figDir,[name '.png']),'Resolution',150);
    savefig(fig,fullfile(figDir,[name '.fig']));
end
