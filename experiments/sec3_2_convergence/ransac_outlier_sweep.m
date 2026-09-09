%% ransac_outlier_sweep.m — 동차다항 RANSAC 의 아웃라이어·수렴 검증 (Sec 3.2)
%
% 세팅: 단일 도형에 occlusion(z-quantile) + 아웃라이어(균일 / 배경평면).
%   동차다항 기저로 적응형 RANSAC(guided sampling + local center refit).
%   반복수는 상수 고정(freezeIter=true). 이론 N은 참고 출력만.
%
% 출력 지표
%   [기하]  cd_occ/cd_vis/cd_all + hd95 + osf (utils/splitSurfaceMetrics, 3.1 과 공용)
%   [분해]  cd_*_oracle = 정답 인라이어로 재적합한 L2, cd_*_direct = RANSAC 없는 L2d
%   [거동]  P/R/F1(점 분류), fp_lev/fn_clump(오분류 기하), fail, gate_margin
%   [자세]  e_par/e_perp(중심 편향의 뷰방향 분해), theta_axis(장축각)
%
% 실행: MATLAB R2025a.  워크스페이스 override: quickMode=true; seedOffset=100; run(...)

run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

% quickMode : 시드 3, maxIter 300 스모크 테스트. 수치 인용 금지.
% devMode   : 축소판. 시드 20 / maxIter 는 full 과 동일.
% computeOracle : L2(oracle) 계산. 셀마다 RANSAC 을 한 번 더 돌려 실행시간이 2배가 된다.
% seedOffset : 0 = tuning fold(λ 선택용), 100 = evaluation fold(최종 수치).
if ~exist('quickMode','var')  || isempty(quickMode),  quickMode  = false; end
if ~exist('devMode','var')    || isempty(devMode),    devMode    = false; end
if ~exist('computeOracle','var') || isempty(computeOracle), computeOracle = ~devMode; end
if ~exist('seedOffset','var') || isempty(seedOffset), seedOffset = 0;     end

%% ---- 공통 파라미터 ----
totalN   = 4000;         % 도형 표면 점 수(가림 전)
noiseSig = 0.007;
baseOcc  = 0.30;         % 기준 가림 비율 (occ 축의 기준셀)
conf     = 0.95;
kMult    = 1.4;          % minimal sample = round(kMult*nT). 4차 21점, 6차 39점
tau      = 1.00;         % inlier threshold (대수적 잔차 |f-1| 기준).
                         % 2026-09-08 τ 스윕에서 τ*≈0.7~1.0(4차) / 1.0~1.5(6차).
                         % 이전 값 0.45 는 ε=0 에서도 참 인라이어를 잘라냈다(recall<1).
failFrac = 0.5;          % cd_all > failFrac*bbox대각 → catastrophic failure
centerFailFrac = 0.15;   % 중심오차 > 15%*bbox대각 → lock-on 실패로 간주
                         % (정상 시드는 ~1%, 평면 ε=0.4 lock-on 실패는 ~30%)
gapMin   = 1.15;         % 장축 유효성: GT 고유값비가 이보다 작으면 자세 평가 불가

% order_axis : RANSAC 기저 차수. 셀마다 기저를 바꿔 4차/6차를 같은 축에서 비교한다.
%              6차는 항 28개라 minimal sample 이 21→39 로 늘고 λ basin 도 다르다.
if quickMode
    nSeed        = 3;
    shapes_axis  = {'hexagon'};
    eps_axis     = [0 0.2];
    lam_axis     = [0 1e-2 1e-1];
    lam_axis6    = [0 1e-1 3e-1];
    order_axis   = [4 6];        % 두 차수 경로를 다 태워야 스모크 테스트가 된다
    maxIter_axis = 300;
    baseMaxIter  = 300;
    otype_axis   = {'uniform','plane'};
    occ_axis     = [0.30 0.50];
    nRepeat      = 4;
elseif devMode
    nSeed        = 20;
    shapes_axis  = {'sphere','cylinder','ellipsoid','hexagon','cube'};
    eps_axis     = [0 0.1 0.2 0.3 0.4];
    lam_axis     = [0 1e-3 1e-2 3e-2 1e-1]; % 4차 basin
    lam_axis6    = [0 3e-3 1e-2 3e-2 1e-1]; % 6차 basin (2026-09-08 격자에서 3e-2~1e-1)
    order_axis   = [4 6];
    maxIter_axis = 1000;
    baseMaxIter  = 1000;
    otype_axis   = {'uniform'};
    occ_axis     = [0.30 0.40 0.50 0.60];
    nRepeat      = 10;
else
    nSeed        = 20;
    shapes_axis  = {'sphere','cylinder','hexagon','cube'};
    eps_axis     = [0 0.1 0.2 0.3 0.4];
    lam_axis     = [0 3e-3 1e-2 3e-2 1e-1 1e0];
    lam_axis6    = [0 1e-1 3e-1 1e0];
    order_axis   = 4;
    maxIter_axis = 1000;   % 예산 실험은 maxIterOverride 로. 축으로 두면 전 블록이 2배가 된다
    baseMaxIter  = 1000;
    otype_axis   = {'uniform','plane'};
    occ_axis     = [0.10 0.20 0.30 0.40 0.50 0.60 0.70];  % 3.1 occ_list 와 동일
    nRepeat      = 10;
end
% baseLam / baseLam6 : 각 차수의 기준 λ. τ=1.0 에서 잰 2D(τ×λ) 격자 기준.
%   λ* 는 occ 에 따라 움직인다. eval fold(τ=1.0, hexagon, ε=.2) argmin:
%     occ .3 → 1e-2 / .4 → 3e-2 / .5 → 3e-2 / .6 → 3e-2
%   관심 대역(occ .4~.6)에 맞춰 3e-2 로 고정. occ .3 에서는 1e-2 대비 22% 손해.
baseEps = 0.2; baseOType = 'uniform';
baseLam = 3e-2;  baseLam6 = 3e-2;
if exist('baseLamOverride','var')  && ~isempty(baseLamOverride),  baseLam  = baseLamOverride;  end
if exist('baseLam6Override','var') && ~isempty(baseLam6Override), baseLam6 = baseLam6Override; end
lam6    = baseLam6;   % 6차 post-hoc refit용 λ
% 축 override: 위 모드 블록이 정한 값을 워크스페이스에서 덮어쓴다.
%   예) baseShape='hexagon_long'; shapesOverride={'hexagon','hexagon_long'}; orderOverride=4;
% epsOverride/occOverride 를 단일 값으로 주면 그 축 블록이 기준셀과 중복돼 사라진다
% (한 축만 길게 보고 싶을 때 나머지를 접는 용도).
if ~exist('baseShape','var')     || isempty(baseShape), baseShape = 'hexagon'; end
if exist('shapesOverride','var') && ~isempty(shapesOverride), shapes_axis = shapesOverride; end
if exist('orderOverride','var')  && ~isempty(orderOverride),  order_axis  = orderOverride;  end
if exist('epsOverride','var')    && ~isempty(epsOverride),    eps_axis    = epsOverride;    end
if exist('occOverride','var')    && ~isempty(occOverride),    occ_axis    = occOverride;    end
% baseOcc override: occ×λ 외의 블록들이 앉는 occ. 관심 대역만 볼 때 같이 옮긴다.
if exist('baseOccOverride','var') && ~isempty(baseOccOverride), baseOcc   = baseOccOverride; end
% maxIter override: maxIter 축은 기준셀 하나만 만들므로, occ·ε 과 교차시키려면
% 예산을 통째로 바꿔 런을 나눠 돌린다.
if exist('maxIterOverride','var') && ~isempty(maxIterOverride)
    maxIter_axis = maxIterOverride; baseMaxIter = maxIterOverride(1);
end
if exist('lam4Override','var')   && ~isempty(lam4Override),   lam_axis    = lam4Override;   end
if exist('lam6Override','var')   && ~isempty(lam6Override),   lam_axis6   = lam6Override;   end
% tau_axis: 인라이어 문턱. 차수마다 잔차 분포의 꼬리가 달라 같은 τ 가 같은 분위수가
% 아니다. 벡터로 주면 축이 되고, 스칼라면 전 셀 공통이다.
if exist('tauOverride','var')    && ~isempty(tauOverride), tau_axis = tauOverride;
else,                                                      tau_axis = tau;        end
baseOrder = min(order_axis);   % 기준셀(수렴곡선·6차 post-hoc refit)이 붙는 차수

% opt_axis : 솔버 동작 ablation. 셀마다 ransacPar 를 갈아끼운다.
%   base    = 현행 (ω 갱신 ON, 국소 재적합을 매 iter)
%   noOmega = ω(점별 샘플링 가중) 갱신 동결 → 유도 샘플링의 기여를 잰다
%   loGate  = updateThresh 0.5 → 국소 재적합을 최고점수의 절반 이상일 때만. 속도 레버
%   both    = 위 둘 동시
% freezeW 는 freezeIter=true 인 한 무의미하다(w 는 maxIter 재계산에만 쓰임).
optPresets = struct( ...
    'base',    struct('freezeOmega',false,'updateThresh',0.0), ...
    'noOmega', struct('freezeOmega',true, 'updateThresh',0.0), ...
    'loGate',  struct('freezeOmega',false,'updateThresh',0.5), ...
    'both',    struct('freezeOmega',true, 'updateThresh',0.5));
if exist('optOverride','var') && ~isempty(optOverride), opt_axis = optOverride;
else,                                                   opt_axis = {'base'}; end

%% ---- 기저 ----
% 차수마다 {Funcs, Grads, nT, Nfit, PhiGrid} 를 한 번 만들어 두고 셀마다 골라 쓴다.
% 6차는 order_axis 에 없어도 post-hoc refit 에 쓰이므로 항상 만든다.
gLin = linspace(-1.9, 1.9, 49); [GX,GY,GZ] = meshgrid(gLin,gLin,gLin);
Gxyz = [GX(:),GY(:),GZ(:)];
Basis = struct('order',{},'Funcs',{},'Grads',{},'nT',{},'Nfit',{},'PhiGrid',{});
for o = unique([order_axis, baseOrder, 6])
    Terms_o = homogeneFischerTerms(o);
    [F_o, G_o, nT_o] = makeFuncsGradsStack(Terms_o);
    Basis(end+1) = struct('order',o,'Funcs',F_o,'Grads',G_o,'nT',nT_o, ...
        'Nfit',round(kMult*nT_o),'PhiGrid',calculateFourthOrder(Gxyz,F_o)); %#ok<SAGROW>
end
basisOf = @(o) Basis(find([Basis.order]==o,1));
B6 = basisOf(6);
nSurf = 12000; hdPct = 95;

%% ---- RANSAC 파라미터 템플릿 ----
ransacPar = struct('maxIter',baseMaxIter,'conf',conf,'thresh',tau, ...
    'minInlierRatio',0.65,'updateThresh',0,'locIters',4, ...
    'momentum',0.66,'damping',0.85,'reg',1e-6,'k',kMult,'lambda',baseLam);
ransacPar.mc  = struct('on',true,'saveVarName','','saveMatFile','','metric','basic');
ransacPar.sim = struct('freezeW',false,'freezeIter',true,'freezeOmega',false);

%% ---- OFAT 셀 구성 (shape,eps,lambda,maxIter,otype,occ,order), 중복 제거 ----
% 전체 블록을 order_axis 로 한 겹 감싼다. λ 는 차수마다 basin 이 달라 따로 준다.
cellList = {};
for ordv = order_axis
for tauv = tau_axis
for mItv = maxIter_axis
for oi2 = 1:numel(opt_axis)
    optv = opt_axis{oi2};
    % lamA = λ×도형 격자, lamB = 기준 λ, lam_occ = occ×λ 격자(OLS + 기준 + 이웃 1)
    % lamN 은 basin 안쪽 이웃. 4차 basin 1e-2~3e-2, 6차 3e-2~1e-1.
    if ordv == 4
        lamA = lam_axis;  lamB = baseLam;  lamN = 1e-2;
    else
        lamA = lam_axis6; lamB = baseLam6; lamN = 1e-1;
    end
    if exist('lamOccOverride','var') && ~isempty(lamOccOverride)
        lam_occ = lamOccOverride; lamEps0 = lamOccOverride;
        if exist('lamEps0Override','var') && ~isempty(lamEps0Override)
            lamEps0 = lamEps0Override;   % ε=0 대조군 λ 를 따로 좁힐 때
        end
    else
        if devMode
            lam_occ = unique([0 lamB lamN]);
        else
            lam_occ = lamA(lamA <= max(lamB,1e-1));
        end
        lamEps0 = unique([0 lamB]);
    end

    % ε 축
    for e = eps_axis
        cellList(end+1,:) = {baseShape, e, lamB, mItv, baseOType, baseOcc, ordv, tauv, optv}; %#ok<*SAGROW>
    end
    % λ × 도형 (leave-one-shape-out 홀드아웃)
    for si = 1:numel(shapes_axis)
        for l = lamA
            cellList(end+1,:) = {shapes_axis{si}, baseEps, l, mItv, baseOType, baseOcc, ordv, tauv, optv};
        end
    end
    % 오염 유형 × ε
    for oi = 1:numel(otype_axis)
        for e = eps_axis(eps_axis>0)
            cellList(end+1,:) = {baseShape, e, lamB, mItv, otype_axis{oi}, baseOcc, ordv, tauv, optv};
        end
    end
    % occ × λ (ε = baseEps)
    for oc = occ_axis
        for l = lam_occ
            cellList(end+1,:) = {baseShape, baseEps, l, mItv, baseOType, oc, ordv, tauv, optv};
        end
    end
    % occ × {OLS, λ*} at ε=0 : 순수 가림 대조군
    for oc = occ_axis
        for l = lamEps0
            cellList(end+1,:) = {baseShape, 0, l, mItv, baseOType, oc, ordv, tauv, optv};
        end
    end
end
end
end
end
key = cellfun(@(s,e,l,m,o,c,d,t,q) sprintf('%s|%.3f|%.4g|%d|%s|%.3f|%d|%.3f|%s',s,e,l,m,o,c,d,t,q), ...
    cellList(:,1),cellList(:,2),cellList(:,3),cellList(:,4),cellList(:,5),cellList(:,6), ...
    cellList(:,7),cellList(:,8),cellList(:,9),'uni',0);
[~,ia] = unique(key,'stable'); cellList = cellList(ia,:);
nCell = size(cellList,1);

%% ---- 출력 폴더 ----
outDir = fullfile(fileparts(mfilename('fullpath')),'..','..','test_ridge');
stamp  = sprintf('sec32_single_%s%s%s', datestr(now,'yyyymmdd_HHMM'), ...
    tern(seedOffset>0,sprintf('_fold%d',seedOffset),''), tern(quickMode,'_quick',''));
figDir = fullfile(outDir, stamp); mkdir(figDir);
fprintf('출력 폴더 → %s\n', figDir);
modeName = tern(quickMode,'quick',tern(devMode,'dev','full'));
fprintf('[%s] 셀 %d개 × 시드 %d회 = 실행 %d회 (oracle %s, order %s, seedOffset=%d)\n', ...
    modeName, nCell, nSeed, nCell*nSeed, tern(computeOracle,'ON','OFF'), ...
    mat2str(order_axis), seedOffset);
% 1행당 실측 ≈ 2.6초 (4차, oracle OFF, maxIter=1000). 비용은 항 수 nT 에 대략 비례.
rowSec = arrayfun(@(i) 2.6*(basisOf(cellList{i,7}).nT/15) ...
    *(1+computeOracle)*(cellList{i,4}/1000), 1:nCell);
fprintf('  예상 소요 ≈ %.0f분\n', sum(rowSec)*nSeed/60);

% 런 설정을 폴더에 남긴다. 나중에 CSV 만 보고 축을 역추적하지 않아도 되게.
writeRunConfig(fullfile(figDir,'run_config.txt'), struct( ...
    'script','ransac_outlier_sweep.m', 'stamp',stamp, 'mode',modeName, ...
    'seedOffset',seedOffset, 'nSeed',nSeed, 'nCell',nCell, 'computeOracle',computeOracle, ...
    'baseShape',baseShape, 'baseEps',baseEps, 'baseOcc',baseOcc, 'baseOType',baseOType, ...
    'baseOrder',baseOrder, 'baseLam',baseLam, 'baseLam6',baseLam6, 'baseMaxIter',baseMaxIter, ...
    'tau_axis',tau_axis, 'order_axis',order_axis, 'lam_axis',lam_axis, 'lam_axis6',lam_axis6, ...
    'eps_axis',eps_axis, 'occ_axis',occ_axis, 'maxIter_axis',maxIter_axis, ...
    'shapes_axis',{shapes_axis}, 'otype_axis',{otype_axis}, 'opt_axis',{opt_axis}, ...
    'totalN',totalN, 'noiseSig',noiseSig, 'kMult',kMult, 'gapMin',gapMin, ...
    'failFrac',failFrac, 'centerFailFrac',centerFailFrac), ransacPar);

% -batch 는 stdout 이 버퍼링되므로 셀마다 이 파일에 한 줄 append 하고 즉시 닫는다.
%   확인:  Get-Content test_ridge\<런폴더>\progress.log -Tail 5 -Wait
logPath = fullfile(figDir,'progress.log');
logLine(logPath,'START %s | cells=%d seeds=%d seedOffset=%d order=%s', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS'), nCell, nSeed, seedOffset, mat2str(order_axis));
tRun = tic;

%% ---- 실행 ----
Rows = {};                       % long-form
convBase = struct('key',{},'curve',{});   % 수렴곡선(기준 관련 셀만)
for ci = 1:nCell
    shp = cellList{ci,1}; epsv = cellList{ci,2};
    lam = cellList{ci,3}; mIt  = cellList{ci,4}; otype = cellList{ci,5};
    occ = cellList{ci,6}; ordv = cellList{ci,7}; tauv = cellList{ci,8};
    optv = cellList{ci,9}; opt = optPresets.(optv);
    Bc  = basisOf(ordv);        % 이 셀이 쓸 기저
    isBaseCell = ordv==baseOrder && strcmp(shp,baseShape) && abs(epsv-baseEps)<1e-9 && ...
                 abs(lam-baseLam)<1e-12 && mIt==baseMaxIter && ...
                 strcmp(otype,baseOType) && abs(occ-baseOcc)<1e-9;
    ransacPar.lambda = lam; ransacPar.maxIter = mIt; ransacPar.thresh = tauv;
    ransacPar.sim.freezeOmega = opt.freezeOmega;
    ransacPar.updateThresh    = opt.updateThresh;
    Ntheory = ceil(log(1-conf)/log(max(realmin,1-(1-epsv)^Bc.Nfit)));  % 참고용

    curveAcc = zeros(mIt,1); curveCnt = 0;
    for sd = 1:nSeed
        rng(sd + seedOffset);
        % --- 도형 생성 + view 회전 (GT=가림/노이즈/아웃라이어 없는 clean 전체) ---
        Pbase = makeShapePoints(shp, totalN);
        Rv = randomRotation();
        PcleanFull = Pbase * Rv.';               % GT surface
        bboxDiag   = norm(max(PcleanFull,[],1) - min(PcleanFull,[],1));
        % --- 가림: z 상위 (1-occ) 만 visible (3.1과 동일 규약) ---
        zThr = prctile(PcleanFull(:,3), 100*occ);
        vis  = PcleanFull(:,3) >= zThr;
        Pin  = PcleanFull(vis,:) + noiseSig*randn(nnz(vis),3);
        nIn  = size(Pin,1);
        % --- 아웃라이어 ---
        nOut = round(epsv/max(1-epsv,1e-9) * nIn);
        switch otype
            case 'uniform'   % visible bbox 1.2x 확장 내부 균일
                if nOut > 0
                    lo = min(Pin,[],1); hi = max(Pin,[],1); c=(lo+hi)/2; hw=1.2*(hi-lo)/2;
                    Pout = c + (2*rand(nOut,3)-1).*hw;
                else
                    Pout = zeros(0,3);
                end
            case 'plane'     % 배경 평면(실측 LiDAR에서 가장 흔한 구조적 오염)
                Pout = makePlaneOutliers(mean(PcleanFull,1), Rv, nOut, zThr, ...
                                         bboxDiag, noiseSig);
            otherwise
                error('Unknown outlier type: %s', otype);
        end
        PPm    = [Pin; Pout];
        gtMask = [true(nIn,1); false(size(Pout,1),1)];
        ransacPar.mc.gtMask = gtMask;

        % --- RANSAC (출력 억제). t_solve = 본 적합 벽시계 시간 ---
        tSolve = tic;
        [~] = evalc(['[ResultDisp,Beta,inMask,Log] = ' ...
            'PoliNavigationSolver3_3_FischerRansac_MC(0,PPm,ordv,Bc.nT,Bc.Funcs,Bc.Grads,ransacPar);']);
        t_solve = toc(tSolve);

        % --- 최종 분류 지표 (점 분류. 표면 F-score 아님) ---
        inMask = logical(inMask);
        TP=sum(gtMask&inMask); FP=sum(~gtMask&inMask); FN=sum(gtMask&~inMask);
        prec=TP/max(1,TP+FP); rec=TP/max(1,TP+FN);
        F1=2*prec*rec/max(1e-12,prec+rec);

        % --- 수렴: running-best postScore, 99% 도달 iter ---
        ps = [Log.postScore]'; ps(isnan(ps)) = -inf;
        rb = cummax(ps); finalBest = rb(end);
        nConv = find(rb >= 0.99*finalBest, 1, 'first'); if isempty(nConv), nConv=mIt; end
        if isBaseCell
            curveAcc = curveAcc + max(rb,0); curveCnt = curveCnt+1;
        end

        % --- τ 게이트 진단: 참 인라이어 잔차의 중앙값이 τ 에 얼마나 근접했나 ---
        % Sobolev 페널티가 l=0 항을 면제하므로 λ→∞ 극한은 β=0 이 아니라 최적 구다.
        % gate_margin < 1 이면 λ 곡선을 정규화 효과로 읽어도 된다.
        rGt = abs(calculateFourthOrder(PPm(gtMask,:) - ResultDisp, Bc.Funcs)*Beta - 1);
        r_med_gt = median(rGt); gate_margin = r_med_gt / tauv;

        % --- L3: 실제 RANSAC 결과의 기하 지표 (3분할) ---
        S3 = splitSurfaceMetrics(Beta, Bc.PhiGrid,GX,GY,GZ, ResultDisp, nSurf, ...
                                 PcleanFull, vis, zThr, hdPct);

        % --- L2 (oracle): 같은 솔버에 정답 인라이어만 넣은 상한 ---
        % 1회 적합으로 대체하면 중심이 가시점 평균에 고정돼 사다리 부호가 뒤집힌다.
        Porc = PPm(gtMask,:);
        if computeOracle
            parOracle = ransacPar; parOracle.mc.gtMask = true(size(Porc,1),1);
            [~] = evalc(['[dispO,betaO,~,~] = PoliNavigationSolver3_3_FischerRansac_MC' ...
                '(0,Porc,ordv,Bc.nT,Bc.Funcs,Bc.Grads,parOracle);']);
            S2 = splitSurfaceMetrics(betaO, Bc.PhiGrid,GX,GY,GZ, dispO, nSurf, ...
                                     PcleanFull, vis, zThr, hdPct);
        else
            S2 = struct('cd_occ',NaN,'cd_all',NaN,'osf',NaN);
        end

        % --- L2d: RANSAC 없이 정답 인라이어를 1회만 적합 (중심=가시점 평균) ---
        % L2 − L2d = RANSAC 의 중심 최적화 몫.
        cO = mean(Porc,1);
        bD = regressionFourthOrder(Porc-cO, Bc.Funcs, lam, ordv, 1);
        Sd = splitSurfaceMetrics(bD, Bc.PhiGrid,GX,GY,GZ, cO, nSurf, ...
                                 PcleanFull, vis, zThr, hdPct);

        % --- 오분류의 기하: 개수(F1)가 못 잡는 것 ---
        fpIdx = ~gtMask & inMask;
        if any(fpIdx)
            dFP = nnDist(PPm(fpIdx,:), PcleanFull);
            fp_lev = median(dFP); fp_lev_p95 = prctile(dFP,95);
        else
            fp_lev = 0; fp_lev_p95 = 0;
        end
        fnIdx = gtMask & ~inMask;
        fn_clump = spreadNN(PPm(fnIdx,:)) / spreadNN(PPm(gtMask,:));

        % --- 중심 편향(뷰방향 분해) + 장축 각도 ---
        eVec  = ResultDisp(:).' - mean(PcleanFull,1);
        e_par = eVec(3) / bboxDiag;              % 가림 방향(월드 z) 성분: 부호 유지
        e_perp= norm(eVec(1:2)) / bboxDiag;      % 직교 성분: 산포

        % --- catastrophic failure: cd_all 문턱 또는 중심 이탈 ---
        % cd_all 만으로는 lock-on 실패를 놓친다(배경평면 ε=0.4 에서 20시드 중 1개만 검출).
        centerErr = norm(eVec) / bboxDiag;
        failFlag  = S3.empty || S3.cd_all > failFrac*bboxDiag || centerErr > centerFailFrac;
        % 자세 게이트를 GT와 추정 양쪽에 건다.
        %  gapGt  ≤ gapMin : 물체가 2차 등방(cube/sphere) → 자세 평가 불가
        %  gapFit ≤ gapMin : 적합 등위면이 블롭이라 주축이 임의 → NaN 으로 둔다
        [axGt, gapGt] = principalAxis(PcleanFull);
        gapFit = NaN; theta_axis = NaN; theta_axis_pca = NaN;
        if ~S3.empty
            [axFit, gapFit] = principalAxis(S3.V);
            if gapGt > gapMin
                axPca = principalAxis(Pin);      % 가시점 raw PCA 베이스라인
                theta_axis_pca = acosd(min(1,abs(dot(axPca,axGt))));
                if gapFit > gapMin
                    theta_axis = acosd(min(1,abs(dot(axFit,axGt))));
                end
            end
        end

        % --- 자세(계수 기반): l=2 블록을 사중극 텐서로 보고 고유분해 ---
        % 등위면을 샘플링하지 않으므로 가려진 쪽으로 뻗은 로브가 주축을 뺏지 않는다.
        % 유효성 게이트는 gap 이 아니라 l2frac(사중극 성분의 비중)으로 건다 —
        % cube/sphere 는 l2frac ~ 0.005 로 정보가 없고, gap 은 통과해버린다.
        [axH, gapH, l2frac] = harmonicAxis(Beta, ordv, Bc.Funcs);
        theta_axis_h = NaN;
        if l2frac > 0.05 && gapGt > gapMin
            theta_axis_h = acosd(min(1,abs(dot(axH,axGt))));
        end

        % --- 6차 post-hoc refit (기준셀만): 최종 inlier에 6차 릿지 ---
        cd_occ6=NaN; cd_all6=NaN;
        if isBaseCell
            Pinl = PPm(inMask,:); cc = mean(Pinl,1);
            b6 = regressionFourthOrder(Pinl-cc, B6.Funcs, lam6, 6, 1);
            S6 = splitSurfaceMetrics(b6, B6.PhiGrid,GX,GY,GZ, cc, nSurf, ...
                                     PcleanFull, vis, zThr, hdPct);
            cd_occ6 = S6.cd_occ; cd_all6 = S6.cd_all;
        end

        Rows(end+1,:) = {shp, epsv, lam, mIt, otype, occ, ordv, tauv, optv, t_solve, sd+seedOffset, nIn, size(Pout,1), Ntheory, ...
            prec, rec, F1, nConv, finalBest, ...
            S3.cd_occ, S3.hd_occ, S3.cd_vis, S3.hd_vis, S3.cd_all, S3.hd_all, S3.osf, ...
            double(S3.occ_empty), ...
            S2.cd_occ, S2.cd_all, S2.osf, Sd.cd_occ, Sd.cd_all, ...
            fp_lev, fp_lev_p95, fn_clump, double(failFlag), ...
            r_med_gt, gate_margin, ...
            e_par, e_perp, theta_axis, theta_axis_pca, gapGt, gapFit, ...
            theta_axis_h, l2frac, gapH, ...
            cd_occ6, cd_all6};
    end
    if curveCnt>0
        convBase(end+1) = struct('key',sprintf('%s e%.2f l%.0e m%d %s o%d',shp,epsv,lam,mIt,otype,ordv), ...
            'curve',curveAcc/curveCnt);
    end
    el = toc(tRun); eta = el/ci*(nCell-ci);
    msg = sprintf('[%2d/%2d] o%d tau=%.2f %-7s %-9s occ=%.2f eps=%.2f lam=%-6.0e mIt=%d %-7s (Ntheory=%d)', ...
        ci,nCell,ordv,tauv,optv,shp,occ,epsv,lam,mIt,otype,Ntheory);
    fprintf('  %s  done\n', msg);
    logLine(logPath,'%s  done | 경과 %s, 남은 예상 %s', msg, hms(el), hms(eta));
end
logLine(logPath,'LOOP DONE %s | 총 %s. 이제 CSV/그림 저장 중...', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS'), hms(toc(tRun)));

T = cell2table(Rows, 'VariableNames', {'shape','eps','lambda','maxIter','otype','occ','order','tau','opt','t_solve','seed', ...
    'n_in','n_out','N_theory','precision','recall','F1','n_conv','bestScore', ...
    'cd_occ','hd_occ','cd_vis','hd_vis','cd_all','hd_all','osf','occ_empty', ...
    'cd_occ_oracle','cd_all_oracle','osf_oracle','cd_occ_direct','cd_all_direct', ...
    'fp_lev','fp_lev_p95','fn_clump','fail','r_med_gt','gate_margin', ...
    'e_par','e_perp','theta_axis','theta_axis_pca','axis_gap_gt','axis_gap_fit', ...
    'theta_axis_h','l2frac','axis_gap_h', ...
    'cd_occ6','cd_all6'});
writetable(T, fullfile(figDir,'ransac_outlier_results.csv'));
save(fullfile(figDir,'ransac_outlier_conv.mat'),'convBase','baseMaxIter');
fprintf('CSV 저장 → %s\n', fullfile(figDir,'ransac_outlier_results.csv'));

%% ---- 반복성(Jaccard): 점군을 고정하고 RANSAC 난수만 바꿔 인라이어 집합 안정성 ----
% 주의: 시드마다 점군 자체가 달라지므로 시드 간 Jaccard 는 정의되지 않는다.
% 여기서는 기준셀 점군 하나를 고정하고 RANSAC 을 nRepeat 회 반복한다.
fprintf('반복성 검사 (기준셀 점군 고정, RANSAC %d회)...\n', nRepeat);
rng(1 + seedOffset);
Pbase = makeShapePoints(baseShape, totalN); Rv = randomRotation();
PcleanFull = Pbase * Rv.';
zThr = prctile(PcleanFull(:,3), 100*baseOcc); vis = PcleanFull(:,3) >= zThr;
Pin  = PcleanFull(vis,:) + noiseSig*randn(nnz(vis),3); nIn = size(Pin,1);
nOut = round(baseEps/(1-baseEps)*nIn);
lo=min(Pin,[],1); hi=max(Pin,[],1); c=(lo+hi)/2; hw=1.2*(hi-lo)/2;
PPm = [Pin; c + (2*rand(nOut,3)-1).*hw];
ransacPar.lambda = baseLam; ransacPar.maxIter = baseMaxIter;
ransacPar.mc.gtMask = [true(nIn,1); false(nOut,1)];
Brep = basisOf(baseOrder);
Masks = false(size(PPm,1), nRepeat);
for rr = 1:nRepeat
    [~] = evalc(['[~,~,imk,~] = ' ...
        'PoliNavigationSolver3_3_FischerRansac_MC(0,PPm,baseOrder,Brep.nT,Brep.Funcs,Brep.Grads,ransacPar);']);
    Masks(:,rr) = logical(imk);
end
jac = [];
for a = 1:nRepeat-1
    for b = a+1:nRepeat
        jac(end+1) = sum(Masks(:,a)&Masks(:,b)) / max(1,sum(Masks(:,a)|Masks(:,b)));
    end
end
Trep = table(median(jac), prctile(jac,25), prctile(jac,75), nRepeat, ...
    'VariableNames',{'jaccard_med','jaccard_q25','jaccard_q75','n_repeat'});
writetable(Trep, fullfile(figDir,'ransac_repeatability.csv'));
fprintf('  Jaccard 중앙값 %.3f (IQR %.3f–%.3f)\n', median(jac), prctile(jac,25), prctile(jac,75));

%% ---- 그림 ----
lamAnchor = [3e-2 1e-1 3e-1];   % 3.1에서 쓰던 lambda 범위 (굵게 표시)
% 기존 그림은 축이 겹치지 않도록 기준 차수만 쓴다. order 4 vs 6 비교는 CSV 로 한다.
makeFigures(T, figDir, struct('shape',baseShape,'eps',baseEps,'lam',baseLam, ...
    'maxIter',baseMaxIter,'otype',baseOType,'occ',baseOcc,'order',baseOrder,'tau',tau_axis(1),'opt',opt_axis{1}));
fprintf('그림 저장 → %s\n', figDir);

%% ============================ local functions ============================
% makeShapePoints / randomRotation / meshSurfacePoints 는 utils/ 로 승격됨.
% 3.1(lcurve_demo)과 같은 코드를 써야 사다리 비교가 성립하므로 여기에 다시 두지 말 것.
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function logLine(path, fmt, varargin)
% 한 줄 append 후 즉시 닫는다(버퍼에 남지 않게). -batch 실행의 진행 확인용.
    fid = fopen(path,'a');
    if fid < 0, return; end
    fprintf(fid, [fmt '\n'], varargin{:});
    fclose(fid);
end

function s = hms(sec)
    if ~isfinite(sec), s='--'; return; end
    h=floor(sec/3600); m=floor(mod(sec,3600)/60);
    if h>0, s=sprintf('%dh%02dm',h,m); else, s=sprintf('%dm%02ds',m,floor(mod(sec,60))); end
end

function Pp = makePlaneOutliers(cObj, Rv, nOut, zThr, extent, sig)
% 배경 평면(벽/바닥) 오염. 뷰프레임 x 방향 법선, 물체 옆에 배치한 뒤
% 타깃과 **동일한 z-quantile 절단**을 적용해 가시성 규약을 일치시킨다.
% 평면이 절단선 아래로 대부분 잘리면 nOut 에 못 미칠 수 있으므로 실제 개수를 반환.
    Pp = zeros(0,3);
    if nOut <= 0, return; end
    nrm = [1 0 0]*Rv.';  e1 = [0 1 0]*Rv.';  e2 = [0 0 1]*Rv.';
    % 거리는 0.5*extent(물체 경계에 스치는 배경)로 둔다. 0.75 에서도 0.5 에서도
    % precision=1.00 이 나오는데, 이는 거리 문제가 아니라 구조적 이유다 — 닫힌
    % 등위면 f=1 은 무한 평면을 흡수할 수 없다. 즉 배경 평면은 균일 아웃라이어보다
    % **쉬운** 오염이며, 이 순서(평면 < 균일 < 2물체)가 §3.2 의 결과 중 하나다.
    cP  = cObj + 0.5*extent*nrm;
    L   = extent;                          % 물체 크기의 ~2배 패치
    tries = 0;
    while size(Pp,1) < nOut && tries < 20
        m = max(4*(nOut - size(Pp,1)), 256);
        Q = cP + (2*rand(m,1)-1)*L.*e1 + (2*rand(m,1)-1)*L.*e2;
        Q = Q + sig*randn(size(Q));
        Pp = [Pp; Q(Q(:,3) >= zThr, :)];  %#ok<AGROW>
        tries = tries + 1;
    end
    if size(Pp,1) > nOut, Pp = Pp(1:nOut,:); end
end

function makeFigures(T, figDir, B)
% Figures for the Sec 3.2 sweep.
% Each figure declares the axis it needs and is skipped when that axis has fewer
% than two levels in T, so a focused run does not emit degenerate plots.
% B: struct of base-cell values (shape, eps, lam, maxIter, otype, occ, order, tau).

    aggMed = @(m,v) median(T.(v)(m & T.fail==0),'omitnan');   % fail 제외
    aggAll = @(m,v) median(T.(v)(m),'omitnan');               % theta 는 전체 시드
    failR  = @(m) 100*mean(T.fail(m));

    epss   = unique(T.eps)';      lams   = unique(T.lambda)';
    occs   = unique(T.occ)';      taus   = unique(T.tau)';
    orders = unique(T.order)';    mits   = unique(T.maxIter)';
    shapes = unique(T.shape);     otypes = unique(T.otype);
    opts   = unique(T.opt);

    % λ=0 은 로그축에 못 올리므로 최소 양수 λ 의 1/3 자리에 놓는다
    lamPos = lams; pos = lams(lams>0);
    if ~isempty(pos), lamPos(lams==0) = min(pos)/3; end
    lamTick = arrayfun(@lamLbl, lams, 'uni', 0);

    %% fig1: inlier classification vs contamination
    if numel(epss) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        nP = 2 + (numel(shapes) >= 2);
        subplot(1,nP,1); hold on; grid on; cols = lines(numel(lams));
        for i = 1:numel(lams)
            y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e,'lambda',lams(i)),'F1'), epss);
            plot(100*epss, y, '-o', 'Color',cols(i,:), 'LineWidth',1.4, 'DisplayName',lamLbl(lams(i)));
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Inlier F1 (median)');
        title('Inlier classification vs contamination');
        legend('Location','southwest','FontSize',7); ylim([0 1]);

        subplot(1,nP,2); hold on; grid on;
        for v = {'precision','recall'}
            y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),v{1}), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',v{1});
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Median');
        title(sprintf('Precision / recall  (\\lambda = %s)', lamLbl(B.lam)));
        legend('Location','southwest','FontSize',7); ylim([0 1]);

        if nP == 3
            subplot(1,3,3); hold on; grid on;
            for si = 1:numel(shapes)
                y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e,'shape',shapes{si}),'F1'), epss);
                plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',shapes{si});
            end
            xlabel('Outlier ratio \epsilon [%]'); ylabel('Inlier F1');
            title('By shape (increasing symmetry order)');
            legend('Location','southwest','FontSize',7); ylim([0 1]);
        end
        saveFig(fig, figDir, 'fig1_classification_vs_eps');
    end

    %% fig2: surface reconstruction error vs contamination
    if numel(epss) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        subplot(1,3,1); hold on; grid on;
        for v = {'cd_occ','cd_vis','cd_all'}
            y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),v{1}), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',strrep(v{1},'_','\_'));
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Chamfer distance');
        title('Chamfer distance, three regions'); legend('Location','northwest','FontSize',7);

        subplot(1,3,2); hold on; grid on;
        for v = {'hd_occ','hd_vis','hd_all'}
            y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),v{1}), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',strrep(v{1},'_','\_'));
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Hausdorff distance (95th pct)');
        title('Hausdorff-95, three regions'); legend('Location','northwest','FontSize',7);

        subplot(1,3,3); hold on; grid on;
        y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),'osf'), epss);
        plot(100*epss, y, '-o', 'LineWidth',1.4);
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Surface fraction in hidden region');
        title('Reconstructed mass behind the occlusion');
        saveFig(fig, figDir, 'fig2_surface_vs_eps');
    end

    %% fig3: loss attribution ladder (needs the oracle refit)
    if numel(epss) >= 2 && any(~isnan(T.cd_all_oracle))
        fig = figure('Visible','off','Position',[50 50 900 400]);
        sfx = {'cd_occ','cd_all'};
        for k = 1:2
            s = sfx{k};
            subplot(1,2,k); hold on; grid on;
            names = {[s '_direct'], [s '_oracle'], s};
            labs  = {'L2d: single ridge, visible centroid', ...
                     'L2: same solver, true inliers', 'L3: RANSAC'};
            for j = 1:3
                y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),names{j}), epss);
                plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',labs{j});
            end
            xlabel('Outlier ratio \epsilon [%]'); ylabel(strrep(s,'_','\_'));
            title(sprintf('Loss attribution: %s', strrep(s,'_','\_')));
            legend('Location','northwest','FontSize',7);
        end
        saveFig(fig, figDir, 'fig3_loss_ladder');
    end

    %% fig4: geometry of the misclassified points
    if numel(epss) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        subplot(1,3,1); hold on; grid on;
        for v = {'fp_lev','fp_lev_p95'}
            y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),v{1}), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',strrep(v{1},'_','\_'));
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Distance from GT surface');
        title('How far false positives lie off the surface');
        legend('Location','northwest','FontSize',7);

        subplot(1,3,2); hold on; grid on;
        y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),'fn_clump'), epss);
        plot(100*epss, y, '-o', 'LineWidth',1.4);
        yline(1,'k--','Uniform spread', 'HandleVisibility','off');
        xlabel('Outlier ratio \epsilon [%]'); ylabel('fn\_clump');
        title('Clustering of discarded inliers  (<1: a whole patch is lost)');

        subplot(1,3,3); hold on; grid on;
        if numel(otypes) >= 2
            for oi = 1:numel(otypes)
                y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e,'otype',otypes{oi}),'cd_all'), epss);
                plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',otypes{oi});
            end
            legend('Location','northwest','FontSize',7);
            title('Contamination type: uniform vs background plane');
            ylabel('Chamfer distance');
        else
            y = arrayfun(@(e) failR(selMask(T,B,'eps',e)), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4);
            title('Catastrophic failure rate'); ylabel('Failure rate [%]');
        end
        xlabel('Outlier ratio \epsilon [%]');
        saveFig(fig, figDir, 'fig4_misclassification');
    end

    %% fig5: pose and centre bias vs contamination
    if numel(epss) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        subplot(1,3,1); hold on; grid on;
        y = arrayfun(@(e) aggMed(selMask(T,B,'eps',e),'e_par'), epss);
        plot(100*epss, y, '-o', 'LineWidth',1.4); yline(0,'k--', 'HandleVisibility','off');
        xlabel('Outlier ratio \epsilon [%]'); ylabel('e_{par} / bbox diagonal');
        title('Centre bias along the viewing direction');

        subplot(1,3,2); hold on; grid on;
        for v = {'e_par','e_perp'}
            y = arrayfun(@(e) abs(aggMed(selMask(T,B,'eps',e),v{1})), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',v{1});
        end
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Normalised distance');
        title('Directed bias vs isotropic scatter');
        legend('Location','northwest','FontSize',7);

        subplot(1,3,3); hold on; grid on;
        pnames = {'theta_axis','theta_axis_pca'};
        plabs  = {'Fitted level set','Raw PCA on visible inliers'};
        for j = 1:2
            y = arrayfun(@(e) aggAll(selMask(T,B,'eps',e),pnames{j}), epss);
            plot(100*epss, y, '-o', 'LineWidth',1.4, 'DisplayName',plabs{j});
        end
        yline(60,'r--','Chance level (60 deg)', 'HandleVisibility','off');
        xlabel('Outlier ratio \epsilon [%]'); ylabel('Major-axis error [deg]');
        title('Major-axis error'); legend('Location','northwest','FontSize',7);
        saveFig(fig, figDir, 'fig5_pose_vs_eps');
    end

    %% fig6: the lambda curve  (primary regularisation figure)
    if numel(lams) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 760]);
        pans = {'cd_all','Chamfer distance (all)',      false; ...
                'cd_occ','Chamfer distance (hidden)',   false; ...
                'theta_axis','Major-axis error [deg]',  true;  ...
                'gate_margin','median |f-1| / \tau',    false};
        for p = 1:4
            subplot(2,3,p); hold on; grid on;
            for c = 1:numel(occs)
                if pans{p,3}
                    y = arrayfun(@(l) aggAll(selMask(T,B,'occ',occs(c),'lambda',l),pans{p,1}), lams);
                else
                    y = arrayfun(@(l) aggMed(selMask(T,B,'occ',occs(c),'lambda',l),pans{p,1}), lams);
                end
                plot(lamPos, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
            end
            if p == 3, yline(60,'r--','Chance', 'HandleVisibility','off'); end
            if p == 4, yline(1,'r--','Gate saturated', 'HandleVisibility','off'); end
            set(gca,'XScale','log','XTick',lamPos,'XTickLabel',lamTick);
            xlabel('\lambda   (leftmost tick = OLS)'); ylabel(pans{p,2});
            title(pans{p,2}); legend('Location','best','FontSize',6);
        end
        subplot(2,3,5); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(l) failR(selMask(T,B,'occ',occs(c),'lambda',l)), lams);
            plot(lamPos, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        set(gca,'XScale','log','XTick',lamPos,'XTickLabel',lamTick);
        xlabel('\lambda'); ylabel('Failure rate [%]'); title('Catastrophic failure rate');
        legend('Location','best','FontSize',6);

        subplot(2,3,6); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(l) aggMed(selMask(T,B,'occ',occs(c),'lambda',l),'axis_gap_fit'), lams);
            plot(lamPos, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        set(gca,'XScale','log','XTick',lamPos,'XTickLabel',lamTick);
        xlabel('\lambda'); ylabel('\lambda_1 / \lambda_2 of fitted surface');
        title('Collapse of the level set towards a sphere');
        legend('Location','best','FontSize',6);
        sgtitle(sprintf('Regularisation sweep  (%s, \\epsilon = %.2f, \\tau = %.2f, order %d)', ...
            B.shape, B.eps, B.tau, B.order));
        saveFig(fig, figDir, 'fig6_lambda_curve');
    end

    %% fig7: does the optimal lambda move with occlusion?
    if numel(lams) >= 3 && numel(occs) >= 2
        fig = figure('Visible','off','Position',[50 50 900 400]);
        subplot(1,2,1); hold on; grid on;
        best = nan(size(occs));
        for c = 1:numel(occs)
            y = arrayfun(@(l) aggMed(selMask(T,B,'occ',occs(c),'lambda',l),'cd_all'), lams);
            plot(lamPos, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
            [mv,i] = min(y);
            if ~isempty(i) && isfinite(mv), best(c) = lamPos(i); end
        end
        set(gca,'XScale','log','XTick',lamPos,'XTickLabel',lamTick);
        xlabel('\lambda'); ylabel('cd\_all'); title('\lambda curve per occlusion level');
        legend('Location','northwest','FontSize',6);

        subplot(1,2,2); hold on; grid on;
        plot(100*occs, best, '-o', 'LineWidth',1.8);
        set(gca,'YScale','log','YTick',lamPos,'YTickLabel',lamTick);
        xlabel('Occlusion [%]'); ylabel('argmin \lambda');
        title({'Optimal \lambda vs occlusion', '(flat = one frozen \lambda is defensible)'});
        saveFig(fig, figDir, 'fig7_lambda_star_vs_occ');
    end

    %% fig8: robustness to occlusion
    if numel(occs) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 760]);
        pans = {'cd_all','Chamfer distance (all)',     false; ...
                'cd_occ','Chamfer distance (hidden)',  false; ...
                'F1','Inlier F1',                      false; ...
                'theta_axis','Major-axis error [deg]', true};
        for p = 1:4
            subplot(2,3,p); hold on; grid on;
            for i = 1:numel(lams)
                if pans{p,3}
                    y = arrayfun(@(c) aggAll(selMask(T,B,'occ',c,'lambda',lams(i)),pans{p,1}), occs);
                else
                    y = arrayfun(@(c) aggMed(selMask(T,B,'occ',c,'lambda',lams(i)),pans{p,1}), occs);
                end
                plot(100*occs, y, '-o', 'LineWidth',1.4, 'DisplayName',lamLbl(lams(i)));
            end
            if p == 4, yline(60,'r--','Chance', 'HandleVisibility','off'); end
            xlabel('Occlusion [%]'); ylabel(pans{p,2}); title(pans{p,2});
            legend('Location','best','FontSize',6);
        end
        subplot(2,3,5); hold on; grid on;
        for i = 1:numel(lams)
            y = arrayfun(@(c) failR(selMask(T,B,'occ',c,'lambda',lams(i))), occs);
            plot(100*occs, y, '-o', 'LineWidth',1.4, 'DisplayName',lamLbl(lams(i)));
        end
        xlabel('Occlusion [%]'); ylabel('Failure rate [%]'); title('Catastrophic failure rate');
        legend('Location','best','FontSize',6);

        subplot(2,3,6); hold on; grid on;
        for i = 1:numel(lams)
            y = arrayfun(@(c) 100*aggAll(selMask(T,B,'occ',c,'lambda',lams(i)),'occ_empty'), occs);
            plot(100*occs, y, '-o', 'LineWidth',1.4, 'DisplayName',lamLbl(lams(i)));
        end
        xlabel('Occlusion [%]'); ylabel('occ\_empty [%]');
        title({'Seeds with no surface in the hidden region', '(there cd\_occ is a bbox penalty, not a distance)'});
        legend('Location','best','FontSize',6);
        sgtitle(sprintf('Occlusion robustness  (%s, \\epsilon = %.2f, \\tau = %.2f, order %d)', ...
            B.shape, B.eps, B.tau, B.order));
        saveFig(fig, figDir, 'fig8_occlusion_sweep');
    end

    %% fig9: does contamination shift or tilt the occlusion curve?
    if numel(occs) >= 2 && any(epss == 0) && any(epss > 0)
        e1 = max(epss); col = lines(numel(lams));
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        vs = {'cd_all','Chamfer distance (all)'; 'theta_axis','Major-axis error [deg]'};
        for p = 1:2
            subplot(1,3,p); hold on; grid on;
            for i = 1:numel(lams)
                yA = arrayfun(@(c) aggAll(selMask(T,B,'occ',c,'lambda',lams(i),'eps',e1),vs{p,1}), occs);
                y0 = arrayfun(@(c) aggAll(selMask(T,B,'occ',c,'lambda',lams(i),'eps',0),vs{p,1}), occs);
                plot(100*occs, yA, '-o', 'Color',col(i,:), 'LineWidth',1.4, ...
                     'DisplayName',sprintf('%s, \\epsilon=%.2f',lamLbl(lams(i)),e1));
                plot(100*occs, y0, '--s', 'Color',col(i,:), 'LineWidth',1.2, ...
                     'DisplayName',sprintf('%s, \\epsilon=0',lamLbl(lams(i))));
            end
            xlabel('Occlusion [%]'); ylabel(vs{p,2});
            title(sprintf('%s  (dashed: no outliers)', vs{p,2}));
            legend('Location','northwest','FontSize',6);
        end
        subplot(1,3,3); hold on; grid on;
        for i = 1:numel(lams)
            yA = arrayfun(@(c) failR(selMask(T,B,'occ',c,'lambda',lams(i),'eps',e1)), occs);
            y0 = arrayfun(@(c) failR(selMask(T,B,'occ',c,'lambda',lams(i),'eps',0)), occs);
            plot(100*occs, yA, '-o', 'Color',col(i,:), 'LineWidth',1.4, ...
                 'DisplayName',sprintf('%s, \\epsilon=%.2f',lamLbl(lams(i)),e1));
            plot(100*occs, y0, '--s', 'Color',col(i,:), 'LineWidth',1.2, ...
                 'DisplayName',sprintf('%s, \\epsilon=0',lamLbl(lams(i))));
        end
        xlabel('Occlusion [%]'); ylabel('Failure rate [%]');
        title('Occlusion alone vs occlusion + contamination');
        legend('Location','northwest','FontSize',6);
        saveFig(fig, figDir, 'fig9_occ_x_eps');
    end

    %% fig10: inlier threshold tau
    if numel(taus) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        subplot(1,3,1); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(t) aggMed(selMask(T,B,'occ',occs(c),'tau',t),'cd_all'), taus);
            plot(taus, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        xlabel('Inlier threshold \tau'); ylabel('cd\_all');
        title('Surface error vs inlier threshold'); legend('Location','best','FontSize',7);

        subplot(1,3,2); hold on; grid on;
        for v = {'precision','recall'}
            y = arrayfun(@(t) aggMed(selMask(T,B,'tau',t),v{1}), taus);
            plot(taus, y, '-o', 'LineWidth',1.6, 'DisplayName',v{1});
        end
        xlabel('Inlier threshold \tau'); ylabel('Median');
        title({'What \tau buys and costs', 'recall < 1 at \epsilon = 0 means true inliers are rejected'});
        legend('Location','best','FontSize',7);

        subplot(1,3,3); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(t) aggAll(selMask(T,B,'occ',occs(c),'tau',t),'theta_axis'), taus);
            plot(taus, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        yline(60,'r--','Chance', 'HandleVisibility','off');
        xlabel('Inlier threshold \tau'); ylabel('Major-axis error [deg]');
        title('Pose error vs inlier threshold'); legend('Location','best','FontSize',7);
        saveFig(fig, figDir, 'fig10_tau_curve');
    end

    %% fig11: basis order comparison (lambda left free: best cell per order)
    if numel(orders) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        vs = {'cd_all','Chamfer distance (all)', false; ...
              'theta_axis','Major-axis error [deg]', true; ...
              'F1','Inlier F1', false};
        for p = 1:3
            subplot(1,3,p); hold on; grid on;
            for o = orders
                y = nan(size(occs));
                for c = 1:numel(occs)
                    vals = arrayfun(@(l) aggMed(selMask(T,B,'occ',occs(c),'order',o,'lambda',l),'cd_all'), lams);
                    [~,bi] = min(vals);
                    if isempty(bi) || ~isfinite(vals(bi)), continue; end
                    m = selMask(T,B,'occ',occs(c),'order',o,'lambda',lams(bi));
                    if vs{p,3}, y(c) = aggAll(m,vs{p,1}); else, y(c) = aggMed(m,vs{p,1}); end
                end
                plot(100*occs, y, '-o', 'LineWidth',1.6, 'DisplayName',sprintf('order %d',o));
            end
            if p == 2, yline(60,'r--','Chance', 'HandleVisibility','off'); end
            xlabel('Occlusion [%]'); ylabel(vs{p,2}); title(vs{p,2});
            legend('Location','best','FontSize',7);
        end
        sgtitle('Basis order comparison, each cell at its own best \lambda');
        saveFig(fig, figDir, 'fig11_order_compare');
    end

    %% fig13: solver-option ablation (speed vs accuracy)
    if numel(opts) >= 2
        fig = figure('Visible','off','Position',[50 50 1300 400]);
        vs = {'t_solve','Solve time per run [s]', false; ...
              'cd_all','Chamfer distance (all)',  false; ...
              'F1','Inlier F1',                   false};
        for p = 1:3
            subplot(1,3,p); hold on; grid on;
            for k = 1:numel(opts)
                y = arrayfun(@(c) aggMed(selMask(T,B,'occ',c,'opt',opts{k}),vs{p,1}), occs);
                plot(100*occs, y, '-o', 'LineWidth',1.6, 'DisplayName',opts{k});
            end
            xlabel('Occlusion [%]'); ylabel(vs{p,2}); title(vs{p,2});
            if p == 1, set(gca,'YScale','log'); end
            legend('Location','best','FontSize',7);
        end
        sgtitle(sprintf('Solver ablation  (%s, \\epsilon = %.2f, \\tau = %.2f, order %d)', ...
            B.shape, B.eps, B.tau, B.order));
        saveFig(fig, figDir, 'fig13_opt_compare');
    end

    %% fig12: iteration budget
    if numel(mits) >= 2
        fig = figure('Visible','off','Position',[50 50 900 400]);
        subplot(1,2,1); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(m) aggMed(selMask(T,B,'occ',occs(c),'maxIter',m),'cd_all'), mits);
            plot(mits, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        set(gca,'XScale','log','XTick',mits);
        xlabel('Iteration budget'); ylabel('cd\_all');
        title('Accuracy vs iteration budget'); legend('Location','best','FontSize',7);

        subplot(1,2,2); hold on; grid on;
        for c = 1:numel(occs)
            y = arrayfun(@(m) prctile(T.n_conv(selMask(T,B,'occ',occs(c),'maxIter',m)),90)/m, mits);
            plot(mits, y, '-o', 'LineWidth',1.4, 'DisplayName',sprintf('occ %.0f%%',100*occs(c)));
        end
        set(gca,'XScale','log','XTick',mits); ylim([0 1.05]);
        xlabel('Iteration budget'); ylabel('n\_conv p90 / budget');
        title({'Score convergence scales with the budget', 'even where accuracy does not'});
        legend('Location','best','FontSize',7);
        saveFig(fig, figDir, 'fig12_maxiter_curve');
    end
end

function m = selMask(T, B, varargin)
% 기준셀에 고정한 마스크. varargin 으로 축 하나씩 값을 바꾼다.
% 값으로 [] 를 주면 그 축은 자유(모두 포함).
    f = struct('shape',B.shape,'eps',B.eps,'lambda',B.lam,'maxIter',B.maxIter, ...
               'otype',B.otype,'occ',B.occ,'order',B.order,'tau',B.tau,'opt',B.opt);
    for i = 1:2:numel(varargin), f.(varargin{i}) = varargin{i+1}; end
    m = true(height(T),1);
    fn = fieldnames(f);
    for i = 1:numel(fn)
        v = f.(fn{i});
        if isempty(v), continue; end
        col = T.(fn{i});
        if iscell(col), m = m & strcmp(col, v);
        else,           m = m & abs(col - v) < 1e-9; end
    end
end

function writeRunConfig(path, S, ransacPar)
% 런 설정을 key = value 로 덤프한다. 사람이 읽고 build_run_index.m 이 파싱한다.
    fid = fopen(path,'w');
    if fid < 0, return; end
    c = onCleanup(@() fclose(fid));
    fprintf(fid, '# run config\n');
    dump(fid, S, '');
    fprintf(fid, '\n# ransacPar\n');
    dump(fid, rmfield(ransacPar,'mc'), '');
    fprintf(fid, '\n# ransacPar.sim\n');
    dump(fid, ransacPar.sim, 'sim.');
end

function dump(fid, S, pre)
    fn = fieldnames(S);
    for i = 1:numel(fn)
        v = S.(fn{i});
        if isstruct(v), continue; end
        if iscell(v),        str = strjoin(cellfun(@char,v,'uni',0), ',');
        elseif ischar(v),    str = v;
        elseif islogical(v), str = mat2str(v);
        elseif isscalar(v),  str = num2str(v,'%g');
        else,                str = strjoin(arrayfun(@(x) num2str(x,'%g'), v, 'uni',0), ',');
        end
        fprintf(fid, '%s%s = %s\n', pre, fn{i}, str);
    end
end

function s = lamLbl(l), if l==0, s='OLS'; else, s=sprintf('%.0e',l); end, end
function saveFig(fig,figDir,name)
    exportgraphics(fig,fullfile(figDir,[name '.png']),'Resolution',150);
    savefig(fig,fullfile(figDir,[name '.fig']));
end
