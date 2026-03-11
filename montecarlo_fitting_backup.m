syms x y z
order  = 4;
TermsC = homogeneFischerTerms(order);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
N = length(TermsC);
%N2 = length(TermsC);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);


range  = [-60 60 -35 40 -30 28];
range_local  = [-10 10 -7 7 -10 10];
faceColor = [0.95, 0.72, 0.5];
faceAlpha = 0.50; meshDen_Data = 80;
warning('off','MATLAB:nearlySingularMatrix')

%%

PP01 = generateRandomPointsOnHexagonPrism(9800)+randn(9800,3)*0.007;
PP02 = generateRandomPointsOnSurface(4000)+randn(4000,3)*0.006;
PP03 = generateRandomPointsOnCylinder(8001)+randn(8001,3)*0.005;
PP04 = (2*randn(200,3)-1)*2;
PP05 = generateRandomPointsOnIcosahedron(4000);
PP_use = PP03;

%%
ransacPar = struct( ...
  'scoreType' ,'similarity',...
  'conf', 0.96, ...
  'thresh', 0.45, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0, ...   % 로컬 최적화 비활성화
  'locIters', 4, ... 
  'momentum', 0.66,...
  'damping', 0.85, ...
  'reg', 1e-6);

[ResultDisp, Beta, inlierMask, Log] = ...
    PoliNavigationSolver3_Fit(PP_use, order, N, Funcs1, Grads1, ransacPar);

% Log: [preScore, postScore, t_raw_ms, t_local_ms]
fprintf('preScore=%.4f, postScore=%.4f, t_raw=%.1f ms, t_local=%.1f ms\n', ...
        Log(1), Log(2), Log(3), Log(4));