%% run_all_shapes.m — lcurve_demo를 5종 도형에 대해 순차 실행 (Sec 3.1)
% 각 도형은 test_ridge/<stamp>_<shape>/ 폴더를 생성한다.
% 실행 후 aggregate_shapes.m 으로 통합 분석.
%
% 배치 실행:
%   & "C:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch ...
%     "run('experiments/sec3_1_basis/run_all_shapes.m')"

shapes_all = {'hexagon', 'cube', 'cylinder', 'sphere', 'ellipsoid'};

for iShape_all = 1:numel(shapes_all)
    % lcurve_demo가 남긴 워크스페이스를 도형 간 격리
    clearvars -except shapes_all iShape_all
    shapeName = shapes_all{iShape_all}; %#ok<NASGU> lcurve_demo가 사용
    quickMode = false;                  %#ok<NASGU>
    fprintf('\n########## [%d/%d] shape = %s ##########\n', ...
        iShape_all, numel(shapes_all), shapes_all{iShape_all});
    run(fullfile(fileparts(mfilename('fullpath')), 'lcurve_demo.m'));
    close all
end
fprintf('\n모든 도형 완료. aggregate_shapes.m 를 실행하세요.\n');
