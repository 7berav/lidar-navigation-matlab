%% 0) 준비
syms x y z
order_Data  = 4;
TermsC_Data = homogeneFischerTerms(order_Data);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
range_Data  = [-45 45 -15 20 -20 25];
faceColor_Data = [0.95, 0.72, 0.5];
faceAlpha_Data = 0.50; meshDen_Data = 100;

%% 1) 데이터 로드
if exist('DATA4.mat','file')
    S_Data = load('DATA4.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA4.mat 파일이 없습니다.');
end


hasFields = all(isfield(Data_Data, {'Beta','Disp'}));
if ~hasFields, error('Data_Data에 Beta/Disp 필드가 없습니다.'); end
valid = arrayfun(@(s) ~isempty(s.Beta) && ~isempty(s.Disp), Data_Data);
idx = find(valid); if isempty(idx), error('유효 항목 없음'); end
outDir = 'func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
% m 파일들이 있는 폴더를 경로에 추가
if ~exist('outDir','var'), outDir = 'func_bank'; end
if exist(outDir,'dir'), addpath(outDir); end
if ~exist(outDir,'dir')
    warning('outDir가 없습니다: %s (이 인덱스는 스킵됩니다)', outDir);
    % 여기서 return 하거나, 루프에서 continue 하세요.
end


%% radial 방향 
figure(12); clf;
ax12 = gca; hold(ax12, 'on');


nG = 25;
cmap = lines(nG);

hLegend = gobjects(0);
for j = 3:3
    iPick = j;
    mpath = sprintf('f_ISS_%d.mat', iPick);
    % 없으면 스킵
    if ~exist(mpath, 'file')
        warning('mat파일을 경로에서 찾지 못했습니다: %s', mpath);
        continue;
    end
     
    c = Data_Data(iPick).Disp(:);       % 중심 벡터 [cx;cy;cz]
    if numel(c) < 3
        error('Disp 길이가 3이 아닙니다.');
    end
    cx = c(1); cy = c(2); cz = c(3);


    S = load(mpath, 'f_shifted_handle');   % 예: S 안에 f_draw(핸들) 또는 Poly/Coeff/Powers 등
    if ~isfield(S, 'f_shifted_handle') || ~isa(S.f_shifted_handle, 'function_handle')
        warning('i=%d: %s 에 f_shifted_handle 함수핸들이 없습니다.', iPick, mpath);
        continue;
    end
    g = S.f_shifted_handle; % alias
    f_rad = @(X,Y,Z) sqrt((X-cx).^2 + (Y-cy).^2 + (Z-cz).^2) .* ...
                 (1 - g(X,Y,Z).^(-1/order_Data));

    fimplicit3(@(X,Y,Z) f_rad(X,Y,Z) - 1, range_Data, ...
        'FaceColor', faceColor_Data, 'EdgeColor','none', ...
        'FaceAlpha', faceAlpha_Data, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) f_rad(X,Y,Z) - 2, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) f_rad(X,Y,Z) - 5, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) f_rad(X,Y,Z) - 10, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
end
axis(ax12, 'equal');
xlim ([-40.5 40.5]); ylim([-20.5 25.5]); zlim([-30.5 30.5]);   % 요청하신 all plot 축 범위
xlabel(ax12, 'X'); ylabel(ax12, 'Y'); zlabel(ax12, 'Z');
title(ax12, 'ISS Stationary Points (from Data\\_Data)');  % 제목만 다르게
view ([1 +0.2 -0.5]);
grid on;
hold(ax12, 'off');
%% gradient based
figure(14); clf;
 gca; hold( 'on');


nG = 25;
cmap = lines(nG);

hLegend = gobjects(0);
for j = 3:4
    iPick = j;
    mpath = sprintf('f_ISS_%d.mat', iPick);
    % 없으면 스킵
    if ~exist(mpath, 'file')
        warning('mat파일을 경로에서 찾지 못했습니다: %s', mpath);
        continue;
    end
     
    c = Data_Data(iPick).Disp(:);       % 중심 벡터 [cx;cy;cz]
    if numel(c) < 3
        error('Disp 길이가 3이 아닙니다.');
    end
    cx = c(1); cy = c(2); cz = c(3);
    Beta = Data_Data(iPick).Beta(:);                % 계수
    d    = Data_Data(iPick).Disp(:);  

    S = load(mpath, 'f_shifted_handle');   % 예: S 안에 f_draw(핸들) 또는 Poly/Coeff/Powers 등
    if ~isfield(S, 'f_shifted_handle') || ~isa(S.f_shifted_handle, 'function_handle')
        warning('i=%d: %s 에 f_shifted_handle 함수핸들이 없습니다.', iPick, mpath);
        continue;
    end
    g0_sym = TermsC_Data * Beta;                     % g0(x,y,z)
    g_sym  = subs(g0_sym, [x,y,z], [x-d(1), y-d(2), z-d(3)]);  % g(x,y,z)
    
    % 그래디언트(심볼릭) → 함수핸들화
    dg_sym = gradient(g_sym, [x y z]);
    g_handle = matlabFunction(g_sym, 'Vars', [x y z]);
     gx_handle = matlabFunction(dg_sym(1), 'Vars', [x y z]);
     gy_handle = matlabFunction(dg_sym(2), 'Vars', [x y z]);
     gz_handle = matlabFunction(dg_sym(3), 'Vars', [x y z]);
    % 1차(법선) 근사 서명거리 핸들
    d_lin_handle = @(X,Y,Z) ( g_handle(X,Y,Z) - 1 ) ./ ...
        sqrt( gx_handle(X,Y,Z).^2 + gy_handle(X,Y,Z).^2 + gz_handle(X,Y,Z).^2 );


    fimplicit3(@(X,Y,Z) d_lin_handle(X,Y,Z) - 0, range_Data, ...
        'FaceColor', faceColor_Data, 'EdgeColor','none', ...
        'FaceAlpha', faceAlpha_Data, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) d_lin_handle(X,Y,Z) - 0.5, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) d_lin_handle(X,Y,Z) - 1, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
    fimplicit3(@(X,Y,Z) d_lin_handle(X,Y,Z) - 1.5, range_Data, ...
        'FaceColor', [0.55, 0.62, 0.9], 'EdgeColor','none', ...
        'FaceAlpha', 0.2, 'MeshDensity', meshDen_Data);
end
axis('equal');
xlim ([-40.5 40.5]); ylim([-20.5 25.5]); zlim([-30.5 30.5]);   % 요청하신 all plot 축 범위
xlabel('X'); ylabel('Y'); zlabel( 'Z');
title( 'ISS Stationary Points (from Data\\_Data)');  % 제목만 다르게
view ([-0.2 +0.5 1]);
grid on;
hold('off');