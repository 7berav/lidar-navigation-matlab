function C = loadSimCloud(matFile)
% sim/raycast_target.py 가 저장한 점군(.mat) 로드
%
% 입력
%   matFile : out_cutting_sim/<run>/complete.mat 또는 s1_<dir>_h<k>.mat
%
% 출력 C (struct)
%   .P          N x 3 점 (노이즈 포함, 타깃 중심 LVLH 좌표 [m])
%   .P_true     N x 3 노이즈 없는 레이 적중점 (평가 기준용)
%   .N          N x 3 삼각형 법선 (센서 쪽)
%   .label      N x 1 uint32 부품 라벨 (1부터, partNames 순서)
%   .range, .incidence(deg), .view_id   N x 1
%   .sensorPos  V x 3 센서 위치
%   .partNames  1 x K cell
%   .meta       struct (R, voxel, lidar, attitude, R_b2l, center_body, stats ...)
%
% body 좌표로 되돌리기:  Pb = C.P * C.meta.R_b2l + C.meta.center_body(:).'
%   (p_l = R_b2l (p_b - c)  →  p_b' = p_l' R_b2l + c')

S = load(matFile);
C.P         = double(S.P);
C.P_true    = double(S.P_true);
C.N         = double(S.N);
C.label     = uint32(S.label(:));
C.range     = double(S.range(:));
C.incidence = double(S.incidence(:));
C.view_id   = uint32(S.view_id(:));
C.sensorPos = double(S.sensor_pos);

pn = S.part_names;
if ischar(pn), pn = cellstr(pn); end
C.partNames = cellfun(@strtrim, pn(:).', 'UniformOutput', false);

C.meta = jsondecode(S.meta_json);      % 중첩 리스트는 행 단위 행렬로 복원됨
if isfield(C.meta, 'R_b2l') && ~isequal(size(C.meta.R_b2l), [3 3])
    error('loadSimCloud: R_b2l 형식 오류');
end
end
