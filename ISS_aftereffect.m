%%데이터 후처리 분리
%함수 늘리고 줄이기 
%있는 다른 점군 beta를 손봐서 그림그리기용 mat 저장. 
%멀쩡한 점군은 그냥 ISS plot manual fix 갈것
syms x y z
order  = 4;
TermsC = homogeneFischerTerms(order);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
range  = [-50 50 -10 15 -30 20];
faceColor = [0.95, 0.72, 0.5];
faceAlpha = 0.50; meshDen = 130;

%% 1) 데이터 로드
outDir = 'func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
if exist(outDir,'dir'), addpath(outDir); end
if ~exist(outDir,'dir')
    warning('outDir가 없습니다: %s (이 인덱스는 스킵됩니다)', outDir);
    % 여기서 return 하거나, 루프에서 continue 하세요.
end

if exist('DATA4.mat','file')
    S_Data = load('DATA4.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA4.mat 파일이 없습니다.');
end
i = 14;

PPm_use = Data_Data(i).P(:);
DispRAN = Data_Data(i).Disp(:);

BetaRAN = Data_Data(2).Beta(:);
%%

theta = pi/3;       % 60도
c = cos(theta);     % 0.5
s = sin(theta);     % 
R  = [1  0   0;0  c  -s;0  s   c]*[0  -1  0;1  0  0;0  0  1];

s_x = 4;   % 길이배수 (2면 x축 길이 2배)
s_y = 2.5;   % y축 스케일 (변경없음)
s_z = 1.8;   % z축 스케일 (변경없음)
S_inv = diag([1/s_x, 1/s_y, 1/s_z]);

Rt = R.';
M = S_inv * Rt;
%% 그림 그려서 테스트
if ~isempty(BetaRAN)
    f1_RAN = TermsC * BetaRAN;
    Funcs1_RAN = matlabFunction(f1_RAN,'Vars',[x y z]);
    f1_rotated = @(X,Y,Z) Funcs1_RAN( ...
    M(1,1).*X + M(1,2).*Y + M(1,3).*Z, ...
    M(2,1).*X + M(2,2).*Y + M(2,3).*Z, ...
    M(3,1).*X + M(3,2).*Y + M(3,3).*Z );


    %PP_inlier = PPm_use(logical(inlierMaskRAN), :);
    %PP_inlier_unbias = PP_inlier - DispRAN;
    Value_inlier_unbias = f1_rotated(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3))-1;
    figure(2)
    scatter3(PPm_use(:,1),PPm_use(:,2),PPm_use(:,3),3,Value_inlier_unbias(:),'filled');
    hold on
    fimplicit3(@(x,y,z) f1_rotated(x,y,z) - 1,[-20 20 -20 20 -16 16],'FaceColor', [0.90, 0.81, 0.53],'EdgeColor','none','FaceAlpha',0.5,'MeshDensity', 250);
    hold off
    colormap(jet);
    colorbar
    caxis ([-0.5 1.0]);
    xlabel ('X (m)')
    ylabel ('Y (m)')
    zlabel ('Z (m)')
    view([1, 1, 1]);
    axis equal
    axis([-40.5 40.5 -20.5 25.5 -30.5 30.5])


end
%%

for i = 14:17

    if i > numel(Data_Data), fprintf('skip %d (index out of range)\n', i); continue; end
    if isempty(Data_Data(i).Disp), fprintf('skip %d (Disp missing)\n', i); continue; end

    Disp = Data_Data(i).Disp(:);           % [dx dy dz]

    orig_x = M(1,1)*(x - Disp(1)) + M(1,2)*(y - Disp(2)) + M(1,3)*(z - Disp(3));
    orig_y = M(2,1)*(x - Disp(1)) + M(2,2)*(y - Disp(2)) + M(2,3)*(z - Disp(3));
    orig_z = M(3,1)*(x - Disp(1)) + M(3,2)*(y - Disp(2)) + M(3,3)*(z - Disp(3));

    f_sym = TermsC * BetaRAN;
    f_rotated_sym = subs(f_sym, [x,y,z], [orig_x, orig_y, orig_z]);
    f_rotated_raw     = matlabFunction(f_rotated_sym, 'Vars', [x y z]);
    f_shifted_handle = @(X,Y,Z) f_rotated_raw(X, Y, Z);

    % func_bank/f_ISS_%d.mat 로 저장 (mpath 변수명 유지)
    mpath = fullfile(outDir, sprintf('f_ISS_%d.mat', i));
    save(mpath, 'f_shifted_handle', '-v7.3');
    
    fname = fullfile(outDir, sprintf('Funcs_ISS_%d.py', i));
    fprintf('Exporting %s ...\n', fname);
    exportFuncToPython(f_rotated_raw,fname, {'x','y','z'});
   
end