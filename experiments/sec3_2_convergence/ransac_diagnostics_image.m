%% ransac_diagnostics_image.m — 3.2 형상·iter진단·성능 시각화
%
% 대표 케이스를 단일 seed 로 재실행하고 전체 Log 를 보존해 그린다.
%   단일  : hexagon, occ=0.30, 균일 아웃라이어 eps=0.20
%   2물체 : hexagon(타깃)+cube(방해물), adjacent, perobj 가림, eps=0.30
% lamShow 의 λ 들을 같은 점군 위에 겹쳐 비교한다.
%
% 그림: A 형상(omega colormap + 등위면, 잔차 colormap 동반)
%       B iter별 진단(preScore/postScore/w/inlierR/||β||²/mseInlier + P/R/F1)
%       C 성능 vs λ (기존 single sweep CSV)
%
% 실행: MATLAB R2025a. override: quickMode=true; run(...)

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
if ~exist('quickMode','var')||isempty(quickMode), quickMode=false; end
% ---- 워크스페이스 override (무엇을 그릴지) ----
%   shapeRep : 도형              (기본 'hexagon')
%   occShow  : occ 목록          (기본 [0.30 0.50 0.70]) — 'occ' 케이스에서 사용
%   epsRep   : 아웃라이어 비율    (기본 0.20)
%   otypeRep : 'uniform'|'plane' (기본 'uniform')
%   lamShow  : 겹쳐 볼 λ 목록     (기본 [0 1e-2 1e-1])
%   seedRep  : 시드              (기본 1)
%   caseSel  : {'occ','single','twoobj'} 중 그릴 것 (기본 전부)
if ~exist('seedRep','var') ||isempty(seedRep),  seedRep = 1;                 end
if ~exist('shapeRep','var')||isempty(shapeRep), shapeRep= 'hexagon';         end
if ~exist('epsRep','var')  ||isempty(epsRep),   epsRep  = 0.20;              end
if ~exist('otypeRep','var')||isempty(otypeRep), otypeRep= 'uniform';         end
if ~exist('lamShow','var') ||isempty(lamShow),  lamShow = [0 1e-2 1e-1];     end
if ~exist('occShow','var') ||isempty(occShow)
    occShow = tern(quickMode,[0.30 0.70],[0.30 0.50 0.70]);
end
if ~exist('caseSel','var') ||isempty(caseSel),  caseSel = {'occ','single','twoobj'}; end
maxIterRep = 300*quickMode + 1000*~quickMode;
lamTag = arrayfun(@(l) tern(l==0,'OLS',sprintf('%.0e',l)), lamShow, 'uni',0);

%% ---- 공통 ----
totalN=4000; noiseSig=0.007; occ=0.30; conf=0.95; kMult=1.4; tau=0.45;
order4=4; TermsC4=homogeneFischerTerms(order4); [Funcs4,Grads4,nT4]=makeFuncsGradsStack(TermsC4);
gLin=linspace(-1.9,1.9,49); [GX,GY,GZ]=meshgrid(gLin,gLin,gLin);
PhiGrid4=calculateFourthOrder([GX(:),GY(:),GZ(:)],Funcs4);

ransacPar=struct('maxIter',maxIterRep,'conf',conf,'thresh',tau,'minInlierRatio',0.65, ...
    'updateThresh',0,'locIters',4,'momentum',0.66,'damping',0.85,'reg',1e-6,'k',kMult,'lambda',0);
ransacPar.mc=struct('on',true,'saveVarName','','saveMatFile','','metric','basic');
ransacPar.sim=struct('freezeW',false,'freezeIter',true,'freezeOmega',false);

outDir=fullfile(fileparts(mfilename('fullpath')),'..','..','test_ridge');
figDir=fullfile(outDir,sprintf('sec32_diag_%s%s',datestr(now,'yyyymmdd_HHMM'),tern(quickMode,'_quick','')));
mkdir(figDir); fprintf('출력 폴더 → %s\n',figDir);

%% ---- 대표 케이스 생성 ----
% 점군 생성 절차는 스윕(ransac_outlier_sweep)과 동일하다.
%   rng(seedRep) -> makeShapePoints -> randomRotation.
% 주의: rng 를 한 번만 걸고 λ 를 순차 실행하므로 난수 스트림이 첫 λ 이후 갈라진다.
%       패널 간 비교는 유효하나 두 번째 이후 λ 의 수치를 CSV 값으로 인용하면 안 된다.
Scase = struct('name',{},'PP',{},'gt',{},'GT',{},'vis',{},'zThr',{});

mkCase = @(nm,PP,gt,GT,vis,zT) struct('name',nm,'PP',PP,'gt',gt,'GT',GT,'vis',vis,'zThr',zT);

% (A) occ 스윕
if any(strcmp(caseSel,'occ'))
    for oc = occShow
        rng(seedRep);
        Pb=makeShapePoints(shapeRep,totalN); Rv=randomRotation(); Pcf=Pb*Rv.';
        zT=prctile(Pcf(:,3),100*oc); vis=Pcf(:,3)>=zT;
        Pin=Pcf(vis,:)+noiseSig*randn(nnz(vis),3); nIn=size(Pin,1);
        nOut=round(epsRep/max(1-epsRep,1e-9)*nIn);
        switch otypeRep
            case 'uniform'
                lo=min(Pin,[],1); hi=max(Pin,[],1); c=(lo+hi)/2; hw=1.2*(hi-lo)/2;
                Pout=c+(2*rand(nOut,3)-1).*hw;
            case 'plane'
                bb=norm(max(Pcf,[],1)-min(Pcf,[],1));
                nrm=[1 0 0]*Rv.'; e1=[0 1 0]*Rv.'; e2=[0 0 1]*Rv.';
                cP=mean(Pcf,1)+0.5*bb*nrm; Pout=zeros(0,3); tr=0;
                while size(Pout,1)<nOut && tr<20
                    m=max(4*(nOut-size(Pout,1)),256);
                    Q=cP+(2*rand(m,1)-1)*bb.*e1+(2*rand(m,1)-1)*bb.*e2+noiseSig*randn(m,3);
                    Pout=[Pout; Q(Q(:,3)>=zT,:)]; tr=tr+1; %#ok<AGROW>
                end
                if size(Pout,1)>nOut, Pout=Pout(1:nOut,:); end
        end
        Scase(end+1)=mkCase(sprintf('occ%02d',round(100*oc)), ...
            [Pin;Pout],[true(nIn,1);false(size(Pout,1),1)],Pcf,vis,zT); %#ok<AGROW>
    end
end

% (B) 단일 기준셀 (occ=0.30)
if any(strcmp(caseSel,'single'))
    rng(seedRep);
    Pb=makeShapePoints(shapeRep,totalN); Rv=randomRotation(); Pcf=Pb*Rv.';
    zT=prctile(Pcf(:,3),100*occ); vis=Pcf(:,3)>=zT;
    Pin=Pcf(vis,:)+noiseSig*randn(nnz(vis),3); nIn=size(Pin,1);
    nOut=round(0.2/0.8*nIn); lo=min(Pin,[],1); hi=max(Pin,[],1); c=(lo+hi)/2; hw=1.2*(hi-lo)/2;
    Pout=c+(2*rand(nOut,3)-1).*hw;
    Scase(end+1)=mkCase('single',[Pin;Pout],[true(nIn,1);false(nOut,1)],Pcf,vis,zT);
end

% (C) 2물체 (hexagon 타깃 + cube 방해물, adjacent, perobj)
if any(strcmp(caseSel,'twoobj'))
    rng(seedRep);
    epsv=0.3; d=1.6; nA=round(totalN*(1-epsv)); nB=round(totalN*epsv);
    Rv2=randomRotation(); A0=makeShapePoints('hexagon',nA)*Rv2.';
    B0=makeShapePoints('cube',nB)*Rv2.'; B0=B0+d*([1 0 0]*Rv2.');
    zTA=prctile(A0(:,3),100*occ); visA=A0(:,3)>=zTA;
    Avis=A0(visA,:); Bvis=B0(B0(:,3)>=prctile(B0(:,3),100*occ),:);
    PP2=[Avis;Bvis]+noiseSig*randn(size(Avis,1)+size(Bvis,1),3);
    Scase(end+1)=mkCase('twoobj',PP2,[true(size(Avis,1),1);false(size(Bvis,1),1)],A0,visA,zTA);
end
fprintf('케이스 %d개: %s\n',numel(Scase),strjoin({Scase.name},', '));

%% ---- 대표 실행 + 형상/시계열 그림 ----
gapMin=1.15;
for ci=1:numel(Scase)
    S=Scase(ci); PP=S.PP; gt=S.gt; GT=S.GT;
    bb=norm(max(GT,[],1)-min(GT,[],1));
    [axGt,gapGt]=principalAxis(GT);
    Runs=struct('lam',{},'Beta',{},'Disp',{},'omega',{},'Log',{},'txt',{});
    for li=1:numel(lamShow)
        ransacPar.lambda=lamShow(li); ransacPar.mc.gtMask=gt;
        [~]=evalc(['[Disp,Beta,inMask,Log,omega] = ' ...
            'PoliNavigationSolver3_3_FischerRansac_MC(0,PP,order4,nT4,Funcs4,Grads4,ransacPar);']);
        % 패널 제목에 그 셀의 지표를 박는다.
        Sm=splitSurfaceMetrics(Beta,PhiGrid4,GX,GY,GZ,Disp,12000,GT,S.vis,S.zThr,95);
        eV=Disp(:).'-mean(GT,1);
        th=NaN;
        if ~Sm.empty
            [axF,gapF]=principalAxis(Sm.V);
            if gapGt>gapMin && gapF>gapMin, th=acosd(min(1,abs(dot(axF,axGt)))); end
        end
        inM=logical(inMask);
        pr=sum(gt&inM)/max(1,sum(inM)); rc=sum(gt&inM)/max(1,sum(gt));
        txt=sprintf(['CD_{occ} %.3f  CD_{all} %.3f\n' ...
                     '|e| %.3f (e_{par} %+.3f)\n' ...
                     '\\theta_{axis} %s   F1 %.3f'], ...
            Sm.cd_occ,Sm.cd_all,norm(eV)/bb,eV(3)/bb, ...
            tern(isnan(th),'n/a (등방)',sprintf('%.0f°',th)), ...
            2*pr*rc/max(1e-12,pr+rc));
        Runs(li)=struct('lam',lamShow(li),'Beta',Beta,'Disp',Disp,'omega',omega,'Log',Log,'txt',txt);
    end

    % ---- A. 형상: omega colormap + 등위면 (1×3) ----
    figShape(Runs, PP, gt, GT, PhiGrid4,GX,GY,GZ, Funcs4, 'omega', figDir, S.name, lamTag);
    figShape(Runs, PP, gt, GT, PhiGrid4,GX,GY,GZ, Funcs4, 'resid', figDir, S.name, lamTag);

    % ---- B. iter별 진단 ----
    figIterTS(Runs, figDir, S.name, lamTag);
    fprintf('  [%s] 형상·시계열 그림 완료\n', S.name);
end

%% ---- C. 성능 vs λ (기존 single sweep CSV) ----
try
    sd=latestDir(outDir,'sec32_single_');
    T=readtable(fullfile(sd,'ransac_outlier_results.csv'));
    lams=unique(T.lambda)'; mi=mode(T.maxIter);
    % otype 축이 있으므로 uniform 으로 한정한다.
    sel=@(l) strcmp(T.shape,'hexagon')&abs(T.eps-0.2)<1e-9&abs(T.lambda-l)<1e-12 ...
             &T.maxIter==mi&strcmp(T.otype,'uniform');
    f1=arrayfun(@(l) median(T.F1(sel(l)),'omitnan'),lams);
    cd=arrayfun(@(l) median(T.cd_all(sel(l)),'omitnan'),lams);
    lp=lams; lp(lp==0)=min(lams(lams>0))/3;
    fig=figure('Visible','off','Position',[60 60 900 400]);
    yyaxis left;  semilogx(lp,f1,'-o','LineWidth',2); ylabel('F1 (median)'); ylim([0 1]);
    yyaxis right; semilogx(lp,cd,'-s','LineWidth',2); ylabel('CD (median)');
    grid on; xlabel('\lambda (0=OLS, 좌측)'); title('성능 vs \lambda (hexagon, \epsilon=0.2) — L-curve 대체');
    xline(1e-2,'k:','\lambda*=1e-2','LineWidth',1.2,'LabelVerticalAlignment','bottom');
    saveFig(fig,figDir,'perf_vs_lambda');
    fprintf('  성능 vs λ 그림 완료\n');
catch ME
    fprintf('  [건너뜀] 성능 vs λ: %s\n', ME.message);
end
fprintf('그림 저장 → %s\n', figDir);

%% ============================ local functions ============================
function s=tern(c,a,b), if c,s=a; else,s=b; end, end
% isoMesh 는 utils/meshSurfacePoints 와 다르다. 표면 점을 샘플링하지 않고
% 렌더링용 삼각 메시(struct)를 돌려준다. 이름이 겹치면 utils 버전을 가리므로 구분한다.
function fv=isoMesh(beta,PhiGrid,GX,GY,GZ,center)
    F=reshape(PhiGrid*beta,size(GX)); fv=isosurface(GX,GY,GZ,F,1);
    if isempty(fv.vertices), fv=struct('vertices',[],'faces',[]); return; end
    fv.vertices=fv.vertices+center;
end
function figShape(Runs,PP,gt,GT,PhiGrid,GX,GY,GZ,Funcs,mode,figDir,name,lamTag)
    nP=numel(Runs);
    fig=figure('Visible','off','Position',[30 30 700*nP 700]);
    for li=1:numel(Runs)
        subplot(1,nP,li); hold on;
        Beta=Runs(li).Beta; Disp=Runs(li).Disp;
        if strcmp(mode,'omega')
            cv=Runs(li).omega; cmap=parula; ttl='omega(점별 가중)';
        else
            cv=Funcs(PP-Disp)*Beta-1; cmap=jet; ttl='잔차 f-1';
        end
        % 아웃라이어(gt=false): 회색으로 작고 옅게
        scatter3(PP(~gt,1),PP(~gt,2),PP(~gt,3),3,[0.6 0.6 0.6],'filled', ...
            'MarkerFaceAlpha',0.35,'HandleVisibility','off');
        % 원래 도형 점(gt=true): colormap
        scatter3(PP(gt,1),PP(gt,2),PP(gt,3),5,cv(gt),'filled');
        if strcmp(mode,'omega'), cl=[min(cv(gt)) max(cv(gt))]; else, cl=[-0.5 0.5]; end
        fv=isoMesh(Beta,PhiGrid,GX,GY,GZ,Disp);
        if ~isempty(fv.vertices)
            patch(fv,'FaceColor',[0.90 0.81 0.53],'EdgeColor','none','FaceAlpha',0.40);
        end
        colormap(gca,cmap); if diff(cl)>0, clim(cl); end; colorbar;
        axis equal; grid on; view(135,20); camlight; lighting gouraud;
        xlabel X; ylabel Y; zlabel Z;
        if Runs(li).lam==0, ls='OLS'; else, ls=sprintf('\\lambda=%.0e',Runs(li).lam); end
        if isfield(Runs,'txt') && ~isempty(Runs(li).txt)
            title({ls, Runs(li).txt},'FontSize',9);   % 수치를 패널에 직접 표기
        else
            title(ls);
        end
    end
    sgtitle(sprintf('%s — %s (회색=아웃라이어/방해물, 컬러=원래 도형) + 회복 등위면',name,ttl));
    saveFig(fig,figDir,sprintf('shape_%s_%s',name,mode));
end
function figIterTS(Runs,figDir,name,lamTag)
    cols=lines(numel(Runs));
    % Fig-a: preScore/postScore/w/inlierR/betaNorm/mseInlier
    fig=figure('Visible','off','Position',[50 50 1200 800]);
    fld={'preScore','postScore','w','inlierR','betaNorm','mseInlier'};
    ylb={'preScore J=S_j','postScore S_{loc}','w (adaptive)','inlierR','||\beta||^2','mseInlier'};
    for pi=1:6
        subplot(3,2,pi); hold on; grid on;
        for li=1:numel(Runs)
            v=[Runs(li).Log.(fld{pi})]'; it=1:numel(v); idx=~isnan(v)&(~strcmp(fld{pi},'betaNorm')|v>0);
            plot(it(idx),v(idx),'.','Color',cols(li,:),'MarkerSize',4,'DisplayName',lamTag{li});
        end
        if strcmp(fld{pi},'betaNorm'), set(gca,'YScale','log'); end
        ylabel(ylb{pi}); if pi>=5, xlabel('iteration'); end
        if pi==1, legend('Location','best','FontSize',7); title(sprintf('%s iter별 점수·w',name)); end
    end
    saveFig(fig,figDir,sprintf('iterts_%s_scores',name));
    % Fig-b: precision/recall/F1 + running-best F1
    fig=figure('Visible','off','Position',[60 60 720 760]);
    fb={'precision','recall','F1'};
    for pi=1:3
        subplot(3,1,pi); hold on; grid on;
        for li=1:numel(Runs)
            v=[Runs(li).Log.(fb{pi})]'; it=1:numel(v); idx=~isnan(v);
            plot(it(idx),v(idx),'.','Color',cols(li,:),'MarkerSize',4,'DisplayName',lamTag{li});
            if strcmp(fb{pi},'F1')  % running-best F1 (postScore 기준 갱신 추적)
                ps=[Runs(li).Log.postScore]'; f1=[Runs(li).Log.F1]';
                rb=nan(size(f1)); best=-inf; bf=NaN;
                for t=1:numel(ps), if ~isnan(ps(t))&&ps(t)>best, best=ps(t); bf=f1(t); end; rb(t)=bf; end
                plot(1:numel(rb),rb,'-','Color',cols(li,:),'LineWidth',1.5,'HandleVisibility','off');
            end
        end
        ylabel(fb{pi}); ylim([0 1]); if pi==3, xlabel('iteration'); end
        if pi==1, legend('Location','southeast','FontSize',7); title(sprintf('%s P/R/F1 (선=running-best F1)',name)); end
    end
    saveFig(fig,figDir,sprintf('iterts_%s_prf',name));
end
function d=latestDir(outRoot,prefix)
    dd=dir(outRoot); nm={dd([dd.isdir]).name};
    hits=nm(startsWith(nm,prefix)&~contains(nm,'_quick')&~contains(nm,'diag'));
    % 진행 중인 런은 폴더만 먼저 생기고 CSV 는 끝날 때 쓰인다 → CSV 있는 것만 후보로.
    hits=sort(hits(cellfun(@(h) isfile(fullfile(outRoot,h,'ransac_outlier_results.csv')),hits)));
    assert(~isempty(hits),'no full-run dir: %s',prefix); d=fullfile(outRoot,hits{end});
end
function saveFig(fig,figDir,name)
    exportgraphics(fig,fullfile(figDir,[name '.png']),'Resolution',150);
    savefig(fig,fullfile(figDir,[name '.fig']));
end
