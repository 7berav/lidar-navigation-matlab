syms x y z


% Fischer term 정의
order_Data  = 4;
TermsC_Data = homogeneFischerTerms(order_Data);   % symbolic 1×N

% 플롯 설정
range_Data  = [-50 50 -10 15 -30 20];
faceColor_Data  = [0.95, 0.72, 0.5];
faceAlpha_Data  = 0.50;
meshDen_Data    = 130;

if exist('DATA3.mat','file')
    S_Data = load('DATA3.mat');
    varNames = fieldnames(S_Data);
    Data_Data = S_Data.(varNames{1});
else
    error('DATA3.mat 파일이 없습니다.');
end

if exist('ISS_stationary.xyz','file')
    P_all_Data = readmatrix('ISS_stationary.xyz','FileType','text');
else
    P_all_Data = [];
end


hasFields = all(isfield(Data_Data, {'Beta','Disp'}));
if ~hasFields
    error('Data_Data에 Beta/Disp 필드가 없습니다. 필드명을 확인하세요: %s', strjoin(fieldnames(Data_Data), ', '));
end


valid = arrayfun(@(s) ~isempty(s.Beta) && ~isempty(s.Disp), Data_Data);
idx = find(valid);
if isempty(idx)
error('유효한 항목이 없습니다.');
end
figure(14); 
if ~isempty(P_all_Data)
    scatter3(P_all_Data(:,1), P_all_Data(:,2), P_all_Data(:,3), 1, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha',.1);
end
hold on
for i = 1:length(idx)
    iPick = idx(i);
   
    hasBeta = ~isempty(Data_Data(iPick).Beta);
    hasDisp = ~isempty(Data_Data(iPick).Disp);
    if  hasBeta && hasDisp
        Beta_Data = Data_Data(iPick).Beta;        % 15x1 double (예시)
        d_Data    = Data_Data(iPick).Disp(:);   % 1x3 row로 정규화
        f_sym_Data = TermsC_Data * Beta_Data;  % symbolic scalar
        f_raw_Data = matlabFunction(f_sym_Data, 'Vars',[x y z]);

        f_shifted_Data = @(X,Y,Z) f_raw_Data(X - d_Data(1), Y - d_Data(2), Z - d_Data(3));
    
    
        fimplicit3(@(X,Y,Z) f_shifted_Data(X,Y,Z) - 1, range_Data, ...
            'FaceColor', faceColor_Data, 'EdgeColor','none', ...
            'FaceAlpha', faceAlpha_Data, 'MeshDensity', meshDen_Data);
        %추출
        %{
        if isfield(Data_Data, 'P') && isfield(Data_Data, 'InlierMask') ...
                && ~isempty(Data_Data(iPick).P) && ~isempty(Data_Data(iPick).InlierMask)
            Pin_Data = Data_Data(iPick).P(logical(Data_Data(iPick).InlierMask), :);
            scatter3(Pin_Data(:,1), Pin_Data(:,2), Pin_Data(:,3), 3, 'r', 'filled');
        end
        %}
       
    else
        warning('i=%d 데이터가 비어있습니다.', iPick);
    end
end
axis equal; %lighting gouraud; camlight headlight;
material dull;
view([0.2 -1 +0.3])
xlabel X; ylabel Y; zlabel Z;
 axis([-40 40 -10 10 -30 20])
title(sprintf('Shifted surface (i=%d)', iPick));
hold off
set(gcf,'Renderer','opengl');       % 투명도는 OpenGL 권장
set(gca,'SortMethod','depth');      % 깊이 기준 정렬
%set(gca,'TwoSidedLighting','on');   % 양면 조명
% 개별 표면(핸들 hSurf가 있을 때)

hSurf = findobj(gca, 'Type', 'Surface');
set(hSurf,'BackFaceLighting','reverselit');  % 뒷면도 빛 받게

%%
for i = 1:21
    iPick = idx(i);
   
    hasBeta = ~isempty(Data_Data(iPick).Beta);
    hasDisp = ~isempty(Data_Data(iPick).Disp);
    if  hasBeta && hasDisp
        Beta_Data = Data_Data(iPick).Beta;        % 15x1 double (예시)
        d_Data    = Data_Data(iPick).Disp(:);   % 1x3 row로 정규화
        f_sym_Data = TermsC_Data * Beta_Data;  % symbolic scalar
        f_shifted_sym = subs(f_sym_Data, [x, y, z], [x - d_Data(1), y - d_Data(2), z - d_Data(3)]);
        f_raw_Data = matlabFunction(f_shifted_sym, 'Vars',[x y z]);

        fname = sprintf('image_RANSAC/Funcs_ISS%d.py', iPick);
        fprintf('Exporting %s ...\n', fname);
        exportFuncToPython(f_raw_Data,fname,{'x','y','z'});
       
    else
        warning('i=%d 데이터가 비어있습니다.', iPick);
    end
end 
