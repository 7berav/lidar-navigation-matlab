function [P_use, info] = prepareISS(opts)
% ISS 실데이터 전처리 (ISS_normcut_4_modified.m 4~57행, run_pipeline_ISS.m 1절과 동일)
%   ISS_stationary.xyz 크롭(P1~P4) + DATA4.mat 세그먼트 5~10 보강 → 반경 thinning
%
% 입력
%   opts : .seed   난수 시드 (thinning 순서)          (기본 1)
%          .rThin  thinning 반경                      (기본 0.15)
%
% 출력
%   P_use : thinning 후 점군 (buildGraph 입력)
%   info  : .N0 (thinning 전 점 수), .rThin
%
% 같은 seed 면 run_pipeline_ISS.m 과 점 집합이 완전히 같다.

if nargin < 1, opts = struct(); end
seed  = getfield_def(opts, 'seed', 1);
rThin = getfield_def(opts, 'rThin', 0.15);
rng(seed);

xyzFile   = 'ISS_stationary.xyz';
data4File = fullfile('resources', 'DATA4.mat');   % resources/는 MATLAB 예약 폴더명이라 addpath 불가

P  = readmatrix(xyzFile, 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21, :);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21, :);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8, :);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 & P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13, :);

S4 = load(data4File);
v4 = fieldnames(S4);
D4 = S4.(v4{1});
Pi = vertcat(D4(5:10).P);

P_use = [P1(1:7:end,:); P2(1:6:end,:); P3(1:11:end,:); P4(1:9:end,:); Pi(1:10:end,:)];

% 반경 thinning: 무작위 순서로 훑으며 반경 r 내 이웃 억제
N0   = size(P_use,1);
nbrs = rangesearch(P_use, P_use, rThin);
ord  = randperm(N0);
keep = false(N0,1); blocked = false(N0,1);
for t = 1:N0
    i = ord(t);
    if blocked(i), continue; end
    keep(i) = true;
    blocked(nbrs{i}) = true;
end
P_use = P_use(keep,:);
info  = struct('N0', N0, 'rThin', rThin);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
