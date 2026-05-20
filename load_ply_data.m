%% load_ply_data.m
% PLY 파일을 읽어 test_ransac / test_ransac_ridge 인풋용 워크스페이스 변수로 로드
%
% 실행 후 워크스페이스에 생성되는 변수:
%   PP_icp_raw      : 20260407_143824_icp_raw.ply       N×3 double
%   PP_obj          : 20260407_143824_object_only.ply   N×3 double
%   PP_d435         : d435_merged_icp_object_only.ply   N×3 double
%
% 각 변수는 정규화 완료 (centroid=0, max_extent 기준 [-1,1] 스케일)
% test 스크립트에서 바로 PPm_use = PP_obj; 처럼 사용 가능

dataDir = 'DATA';

plyFiles = { ...
    '20260407_143824_icp_raw.ply',          'PP_icp_raw'; ...
    '20260407_143824_object_only.ply',      'PP_obj';     ...
    'd435_merged_icp_object_only.ply',      'PP_d435';    ...
};

for fi = 1:size(plyFiles, 1)
    plyPath = fullfile(dataDir, plyFiles{fi, 1});
    varName = plyFiles{fi, 2};

    if ~isfile(plyPath)
        warning('파일 없음: %s → 건너뜀', plyFiles{fi,1});
        continue;
    end

    fprintf('읽는 중: %s\n', plyFiles{fi,1});

    pc = pcread(plyPath);
    PP = double(pc.Location);   % N×3

    % NaN/Inf 제거
    PP = PP(all(isfinite(PP),2), :);

    % 정규화: centroid=0, max_extent 기준 스케일
    c  = mean(PP, 1);
    PP = PP - c;
    s  = max(abs(PP(:)));
    if s > eps,  PP = PP / s;  end

    % 워크스페이스에 저장
    assignin('base', varName, PP);

    fprintf('  → %s  (%d pts)  X:[%.3f %.3f]  Y:[%.3f %.3f]  Z:[%.3f %.3f]\n', ...
        varName, size(PP,1), ...
        min(PP(:,1)), max(PP(:,1)), ...
        min(PP(:,2)), max(PP(:,2)), ...
        min(PP(:,3)), max(PP(:,3)));
end

fprintf('\n완료. 사용 예시:\n');
fprintf('  PPm_use = PP_obj;\n');
fprintf('  gtMask  = true(size(PPm_use,1),1);\n');

%% ---- 시각화 ----
varNames = plyFiles(:, 2);
n = numel(varNames);

figure('Name', 'PLY Point Clouds', 'Position', [100 100 400*n 500]);
for fi = 1:n
    if ~evalin('base', sprintf('exist(''%s'',''var'')', varNames{fi}))
        continue;
    end
    PP = evalin('base', varNames{fi});

    subplot(1, n, fi);
    scatter3(PP(:,1), PP(:,2), PP(:,3), 1, 'filled');
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(fi);
    view(3);
end


