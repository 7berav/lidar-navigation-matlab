%% 0) 준비
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
%% Data.mat 불러오기
if ~exist('image_cut/DATA7.mat','file')
    error('DATA7.mat 파일이 없습니다. 스펙트럴 분할 결과를 먼저 저장하십시오.');
end
S = load('image_cut/DATA7.mat');           % 기대 변수: Data6, Kmix
if ~isfield(S,'Data')
    error('DATA7.mat 안에 Data 변수가 없습니다.');
end
Data = S.Data;
if isfield(S,'Kmix'), Kmix = S.Kmix; 
else, Kmix = numel(Data); 
end


if exist('ISS_stationary.xyz','file')
    P_all_Data = readmatrix('ISS_stationary.xyz','FileType','text');
else
    P_all_Data = [];
end
outDir = 'image_cut/func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
if exist(outDir,'dir'), addpath(outDir); end
%% Figure 2: 군집별 점 그리기
figure(1); clf; 
scatter3(P_all_Data(:,1), P_all_Data(:,2), P_all_Data(:,3), 2, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.10, 'MarkerEdgeAlpha', 0.10, 'DisplayName', 'P\_all\_Data (bg)');

hold on
cmap = lines(Kmix);
for i = 18:23
    if isempty(Data(i).P), continue; end
    G = Data(i).P;
    scatter3(G(:,1), G(:,2), G(:,3), 4, cmap(i,:), 'filled');
end
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Clusters from Data');
hold off
 
axis('equal');
axis([-40.5 40.5 -20.5 25.5 -30.5 30.5]);   % 요청하신 all plot 축 범위
title( 'ISS Station Points'); 
view([0.2 -1 -0.1])

%% 공통 RANSAC 파라미터
ransacPar = struct('maxIter',6600,'conf',0.96,'thresh',0.400, ...
                   'minInlierRatio',0.65,'updateThresh',0.55, ...
                   'locIters',4,'momentum', 0.66,'damping',0.85,'reg',1e-6);

%figure(3); clf
%scatter3(P_all_Data(:,1), P_all_Data(:,2), P_all_Data(:,3), 1, 'k', 'filled'); hold on
    
for k = 16:Kmix
    Pk = Data(k).P;
    if isempty(Pk), continue; end

    % RANSAC 실행
    PPm_use = Pk;
    [DispRAN, BetaRAN, inlierMaskRAN] = PoliNavigationSolver3_FischerRansac( ...
        0, PPm_use, order, length(TermsC), Funcs1, Grads1, ransacPar);

    eval(sprintf('DispRAN%d = DispRAN;', k));
    eval(sprintf('BetaRAN%d = BetaRAN;', k));
    eval(sprintf('inlierMaskRAN%d = inlierMaskRAN;', k));
    
    Data(k).Disp        = DispRAN;              % 1x3 or 3x1
    Data(k).Beta        = BetaRAN(:);
    Data(k).InlierMask  = logical(inlierMaskRAN(:));

    if ~isempty(BetaRAN)
        eval(sprintf('PPm_use_unbias%d = PPm_use - DispRAN;', k));
        
    end
end



%% graph for any ransac

%target setting (whether if 4 or 5)
k_draw = 14;
PPm_use = Data(k_draw).P;
%BetaRAN = BetaRAN7;
%DispRAN = DispRAN7;
%inlierMaskRAN = inlierMaskRAN7;
%PPm_use_unbias= PPm_use_unbias7;
eval(sprintf('DispRAN = DispRAN%d;', k_draw));
eval(sprintf('BetaRAN = BetaRAN%d;', k_draw));
eval(sprintf('inlierMaskRAN = inlierMaskRAN%d;', k_draw));
eval(sprintf('PPm_use_unbias = PPm_use_unbias%d;', k_draw));

if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN);

    PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = Funcs1_RAN(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3))-1;
    figure(2)
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(f1_RAN-1,range_local,'FaceColor', [0.90, 0.81, 0.53], ...
        'EdgeColor','none','FaceAlpha',0.5,'MeshDensity',meshDen_Data);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 0.5]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1 -1 1])
    axis equal
    axis([-10.5 10.5 -10.5 10.5 -10.5 10.5]);
    
    %axis(range_local)

    
    figure(3)
    scatter3(PPm_use_unbias(:,1),PPm_use_unbias(:,2),PPm_use_unbias(:,3),2,[0.5 0.5 0.5],'filled');
    hold on
    scatter3(PP_inlier_unbias(:,1),PP_inlier_unbias(:,2),PP_inlier_unbias(:,3),3,'b','filled');
    
    fimplicit3(f1_RAN-1,range_local,'FaceColor', [0.95, 0.82, 0.5], ...
        'EdgeColor','none','FaceAlpha',0.5, 'MeshDensity', meshDen_Data);
    hold off
    colormap(jet);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1 -1 1])
    axis equal
    axis([-10.5 10.5 -10.5 10.5 -10.5 10.5]);
    %axis equal
    %axis(range_local)

end
%% data sace
matPath = 'image_cut/DATA7.mat';
if ~exist(matPath,'file')
    error('DATA7.mat 파일이 없습니다. 먼저 생성하세요.');
end
M = matfile(matPath, 'Writable', true);

% Data 전체 교체
M.Data = Data;                 % 구조체 배열 Data를 통째로 갱신
M.Kmix = numel(Data);   



%% function save

%{
for k = 11:Kmix    
    iPick = k;
    Beta_Data = Data(iPick).Beta(:);
    d_Data    = Data(iPick).Disp(:);
  % 원점 기준 다항식
    f_sym_Data = TermsC * Beta_Data;

    % 평행이동만 반영
    f_shifted_sym = subs(f_sym_Data, [x, y, z], [x - d_Data(1), y - d_Data(2), z - d_Data(3)]);

    % 함수핸들화(그림/내보내기 공통 베이스)
    f_shifted_handle = matlabFunction(f_shifted_sym, 'Vars', [x y z]);   % @(x,y,z)
    % mat 저장 
    %expr = char(f_shifted_sym);     % .* ./ .^ 유지
    save(fullfile(outDir, sprintf('f_ISS_%d.mat', iPick)), 'f_shifted_handle', '-v7');

    % i별 파일로 저장(요청하신 exportFuncToPython 그대로 사용)
    fname_py = fullfile(outDir, sprintf('Funcs_ISS_%d.py', iPick));
    fprintf('Exporting (shift-only) %s ...\n', fname_py);
    exportFuncToPython(f_shifted_handle, fname_py, {'x','y','z'});
end

%}

%{
Data3(25) = struct('P',[],'Disp',[],'Beta',[],'InlierMask',[]);

for i = 1:25
    vB = sprintf('BetaRAN%d', i);
    if ~exist(vB,'var')
        fprintf('skip %d (Beta missing)\n', i); continue;
    end
    Beta = eval(vB);
    if isempty(Beta)
        fprintf('skip %d (Beta empty)\n', i); continue;
    end

    vD = sprintf('DispRAN%d', i);
    if ~exist(vD,'var'), fprintf('skip %d (Disp missing)\n', i); continue; end
    Disp = eval(vD);  Disp = Disp(:);  % [dx dy dz] 행벡터 보장
    vP = sprintf('P%d', i);
    if exist(vP,'var'), P = eval(vP); else, P = []; end
    vM = sprintf('inlierMaskRAN%d', i);
    if exist(vM,'var'), M = eval(vM); else, M = []; end

    Data3(i).P          = P;
    Data3(i).Disp       = Disp;
    Data3(i).Beta       = Beta;
    Data3(i).InlierMask = M;
    f_raw_sym = TermsC * Beta;
    f_raw     = matlabFunction(f_raw_sym, 'Vars', [x y z]);
    f_shifted = @(X,Y,Z) f_raw(X - Disp(1), Y - Disp(2), Z - Disp(3));
    
    fname = sprintf('%s/Funcs_ISS%d.py', outDir, i);
    fprintf('Exporting %s ...\n', fname);
    exportFuncToPython({f_shifted},fname, 1, {}, {'x','y','z'});
   
end

%}
