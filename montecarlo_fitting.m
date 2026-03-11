syms x y z
order  = 4;
TermsC = homogeneFischerTerms(order);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
nT = length(TermsC);
%N2 = length(TermsC);
[Funcs1, Grads1, ~]=makeFuncsGradsStack(TermsC);

rot_each   = true; 
range  = [-60 60 -35 40 -30 28];
range_local  = [-10 10 -7 7 -10 10];
faceColor = [0.95, 0.72, 0.5];
faceAlpha = 0.50; meshDen_Data = 80;
warning('off','MATLAB:nearlySingularMatrix')

%%

PP01 = generateRandomPointsOnHexagonPrism(10000);
PP02 = generateRandomPointsOnSurface(4000)+randn(4000,3)*0.006;
PP03 = generateRandomPointsOnCylinder(10000);
PP04 = (2*randn(200,3)-1)*2;
PP05 = generateRandomPointsOnIcosahedron(4000);
PP_use = PP03;

if rot_each
    % 간단 랜덤 직교행렬 (QR 기반)
    A = randn(3); [Q,~] = qr(A);
    PP_use = (Q * PP_use.').';   % Nx3
end
%%
fitPar = struct( ...
  'scoreType' ,'similarity',...
  'conf', 0.96, ...
  'thresh', 0.45, ...
  'minInlierRatio', 0.65, ...
  'updateThresh', 0, ...   % 로컬 최적화 비활성화
  'locIters', 10, ... 
  'momentum', 0.66,...
  'damping', 0.85, ...
  'reg', 1e-6);

[ResultDisp, Beta, inlierMask, Log] = ...
    PoliNavigationSolver3_Fit(PP_use, order, nT, Funcs1, Grads1, fitPar);

% Log: [preScore, postScore, t_raw_ms, t_local_ms]
fprintf('preScore=%.4f, postScore=%.4f, t_raw=%.1f ms, t_local=%.1f ms\n', ...
        Log(1), Log(2), Log(3), Log(4));

%%
% === 5) 결과 버퍼 ===
Res = [];   % 테이블로 모을 예정
row = 0;
%%

sigma_list = [ 0.05 ];     % 가우시안 노이즈 표준편차
q_list     = [0.70 0.5 0.0 0.1 0.30];  
%%
% === 6) 실험 루프 ===
for s = sigma_list
  % 노이즈 추가
  Pnoisy = PP_use + randn(size(PP_use))*s;

  for q = q_list
    P = Pnoisy;



    % 가시율: z-축 프로젝션 상위 q 분위수만 유지
    zproj = P(:,3);
    tau   = quantile(zproj, q);    % 또는 prctile(zproj, q*100)
    keep  = zproj >= tau;
    P_use = P(keep, :);

    if size(P_use,1) < nT*2
        fprintf('[skip] sigma=%.3f, q=%.2f → 표본 너무 적음(N=%d)\n', s, q, size(P_use,1));
        continue;
    end

    % 단일 피팅 실행 (회귀1회 + 로컬1회)
    [Disp, Beta, inlierMask, Log] = ...
        PoliNavigationSolver3_Fit(P_use, order, nT, Funcs1, Grads1, fitPar);
    
    % Log = [preScore, postScore, raw_ms, local_ms]
    % 결과 집계
    row = row + 1;
    Res(row).sigma     = s;

    Res(row).q         = q;
    Res(row).iter      = fitPar.locIters;    
    Res(row).N         = size(P_use,1);
    Res(row).preScore  = Log(1);
    Res(row).postScore = Log(2);
    Res(row).t_raw_ms  = Log(3);
    Res(row).t_loc_ms  = Log(4);
    Res(row).inlierR   = mean(inlierMask);   % 참고치
    
    fprintf('σ=%.3f, q=%.2f, N=%5d | pre=%.4f → post=%.4f | raw=%.1fms, loc=%.1fms\n', ...
        s, q, Res(row).N, Res(row).preScore, Res(row).postScore, Res(row).t_raw_ms, Res(row).t_loc_ms);
  
    % 입력: P_use (Nx3), Beta (nT×1), ResultDisp(1×3), Funcs, Grads
    P_result = P - Disp;                          % 모델 좌표계로 평행이동
    phi   = Funcs1([P_result(:,1), P_result(:,2), P_result(:,3)]);   % N×nT
    Fval  = phi * Beta - 1;                              % N×1, 암시적 잔차
    
    Gstack = Grads1(P_result);                               % (3N)×nT
    GradVal = Gstack * Beta;                             % (3N)×1
    Nn = size(P_result,1);
    gn = sqrt( max(1e-12, ...
          GradVal(1:Nn).^2 + GradVal(Nn+1:2*Nn).^2 + GradVal(2*Nn+1:3*Nn).^2) );
    
    d_sdf = Fval ./ gn;                % 부호 거리
    %d_sdf = abs(Fval) ./ gn;                % 절대 거리

    SDF_RMS = sqrt(mean(d_sdf.^2));
    SDF_P95 = prctile(d_sdf, 95);
    Res(row).Dispx     = Disp(1);
    Res(row).Dispy     = Disp(2);
    Res(row).Dispz     = Disp(3);
    Res(row).SDF_RMS   = SDF_RMS;
    Res(row).SDF_P95   = SDF_P95;
  end
end

% === 7) 표로 정리 ===
T = struct2table(Res);
%%
d = d_sdf;

figure(21); clf;
histogram(d, 'NumBins', 80);    % 필요하면 'Normalization','pdf'
grid on; box on;
xlabel('d\_sdf'); ylabel('count');
title('Histogram of d\_sdf');


%%
PPm_use_unbias =  PP_use - Disp;

if ~isempty(Beta)
    f1_RAN = TermsC * Beta;
    Funcs1_RAN = matlabFunction(f1_RAN);
    
    PP_inlier = P_use;
    %PP_inlier = P_use(logical(inlierMask), :);
    PP_inlier_unbias = PP_inlier - Disp;
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
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

    
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
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])

end