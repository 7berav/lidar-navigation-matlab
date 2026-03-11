%% 0) 준비
syms x y z
order_Data  = 4;
TermsC_Data = homogeneFischerTerms(order_Data);   % (Fisher 기반 생성 가능. 비-Fisher면 직접식 대입)
range_Data  = [-60 60 -35 40 -30 28];
faceColor_Data = [0.95, 0.72, 0.5];
faceAlpha_Data = 0.50; meshDen_Data = 200;

%% 1) 데이터 로드
if exist('image_cut/DATA7.mat','file')
    S_Data = load('image_cut/DATA7.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
    K = S_Data.Kmix;
else
    error('DATA7.mat 파일이 없습니다.');
end

if exist('ISS_stationary.xyz','file')
    P_all_Data = readmatrix('ISS_stationary.xyz','FileType','text');
else
    P_all_Data = [];
end

hasFields = all(isfield(Data_Data, {'Beta','Disp'}));
if ~hasFields, error('Data_Data에 Beta/Disp 필드가 없습니다.'); end
valid = arrayfun(@(s) ~isempty(s.Beta) && ~isempty(s.Disp), Data_Data);
idx = find(valid); if isempty(idx), error('유효 항목 없음'); end


%디렉토리에 이 폴더를 만들어서 여기에 저장하면 좋은데, 절대 좌표 아니라 addpath로 하는거라 어디에 있어도 다 동작함. 저장은 
%outDir = 'func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
outDir = 'image_cut/func_bank'; if ~exist(outDir,'dir'), mkdir(outDir); end
if exist(outDir,'dir'), addpath(outDir); end
%% 점만 그리기?
figure(12); clf;
ax12 = gca; hold(ax12, 'on');

   %scatter3(ax12, P_all_Data(:,1), P_all_Data(:,2), P_all_Data(:,3), 2, 'k', 'filled', ...
      %  'MarkerFaceAlpha', 0.10, 'MarkerEdgeAlpha', 0.10, 'DisplayName', 'P\_all\_Data (bg)');



nG = Kmix;
cmap = lines(nG);

hLegend = gobjects(0);
for j = 1:nG
    iPick = j;

    % 점군 읽기: Data_Data(iPick).P (없으면 스킵)
    if ~isfield(Data_Data, 'P') || isempty(Data_Data(iPick).P)
        continue;
    end
    G = Data_Data(iPick).P;
    if size(G,2) < 3 || isempty(G)
        continue;
    end

    % (선택) 인라이어 마스크가 있으면 적용
    if isfield(Data_Data, 'InlierMask') && ~isempty(Data_Data(iPick).InlierMask)
        m = logical(Data_Data(iPick).InlierMask);
        if numel(m) == size(G,1)
            G = G(m, :);
        end
    end

    % 색상 선택 및 그리기
    c = cmap(j, :);
    h = scatter3(ax12, G(:,1), G(:,2), G(:,3), 6, c, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85, ...
        'DisplayName', sprintf('Data\\_Data(%d)', iPick));
    hLegend(end+1) = h; %#ok<SAGROW>
end
axis(ax12, 'equal');
axis([-60 60 -40 40 -30 30])  % 요청하신 all plot 축 범위
xlabel(ax12, 'X'); ylabel(ax12, 'Y'); zlabel(ax12, 'Z');
view([0.2 -1 -0.1])
title(ax12, 'ISS Stationary Points (from Data\\_Data)');  % 제목만 다르게
hold(ax12, 'off');
%% 존재하는 케이스만 함수 mat 저장 -> 한번만 돌려서 mat 있으면 안해도 됨
%여기서 오류떠서 손수해야되는 애들은 aftereffect 코드 찾아보도록
%{
for k = 1:numel(idx)
    iPick = idx(k);
    Beta_Data = Data_Data(iPick).Beta(:);
    d_Data    = Data_Data(iPick).Disp(:);
  % 원점 기준 다항식
    f_sym_Data = TermsC_Data * Beta_Data;

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

%% loop 2 graph
%함수 .mat 읽어와서 그리기. mat은 
figure(14); %clf;
%if ~isempty(P_all_Data)
%    scatter3(P_all_Data(:,1), P_all_Data(:,2), P_all_Data(:,3), 1, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',.1);
%end
hold on

% m 파일들이 있는 폴더를 경로에 추가
if ~exist('outDir','var'), outDir = 'func_bank'; end
if exist(outDir,'dir'), addpath(outDir); end
if ~exist(outDir,'dir')
    warning('outDir가 없습니다: %s (이 인덱스는 스킵됩니다)', outDir);
    % 여기서 return 하거나, 루프에서 continue 하세요.
end

for iPick = 1:25
    mpath = sprintf('f_ISS_%d.mat', iPick);
    % 없으면 스킵
    if ~exist(mpath, 'file')
        warning('mat파일을 경로에서 찾지 못했습니다: %s', mpath);
        continue;
    end
     
    S = load(mpath, 'f_shifted_handle');   % 예: S 안에 f_draw(핸들) 또는 Poly/Coeff/Powers 등
    % 존재 확인
    if ~isfield(S, 'f_shifted_handle') || ~isa(S.f_shifted_handle, 'function_handle')
        warning('i=%d: %s 에 f_shifted_handle 함수핸들이 없습니다.', iPick, mpath);
        continue;
    end

    f_draw =S.f_shifted_handle;

    

    % 그리기 (등표면 f=1)
    fimplicit3(@(X,Y,Z) f_draw(X,Y,Z) - 1, range_Data, ...
        'FaceColor', faceColor_Data, 'EdgeColor','none', ...
        'FaceAlpha', faceAlpha_Data, 'MeshDensity', meshDen_Data);

    % (선택) 인라이어 포인트 표시 — 필요 시 주석 해제
    % if isfield(Data_Data, 'P') && isfield(Data_Data, 'InlierMask') ...
    %        && ~isempty(Data_Data(iPick).P) && ~isempty(Data_Data(iPick).InlierMask)
    %     Pin_Data = Data_Data(iPick).P(logical(Data_Data(iPick).InlierMask), :);
    %     scatter3(Pin_Data(:,1), Pin_Data(:,2), Pin_Data(:,3), 3, 'r', 'filled');
    % end
end
axis equal; lighting gouraud; camlight headlight;
material dull;
view([1 -1 -1])
xlabel X; ylabel Y; zlabel Z;
 axis([-60 60 -40 40 -30 30])
title(sprintf('ISS surface'));
hold off
set(gcf,'Renderer','opengl');       % 투명도는 OpenGL 권장
set(gca,'SortMethod','depth');      % 깊이 기준 정렬
%set(gca,'TwoSidedLighting','on');   % 양면 조명
% 개별 표면(핸들 hSurf가 있을 때)

hSurf = findobj(gca, 'Type', 'Surface');
set(hSurf,'BackFaceLighting','reverselit');  % 뒷면도 빛 받게

