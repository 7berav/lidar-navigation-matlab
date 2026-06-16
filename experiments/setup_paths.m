% setup_paths.m  Add project paths for paper/2026taes experiments.
%
% Invoke from any experiment script (run as a script, not a function):
%   run(fullfile(fileparts(mfilename('fullpath')), '../../experiments/setup_paths.m'))
%
% Adds the repository root, utils/, and experiments/ (recursively) to the
% MATLAB path for the current session. Idempotent.

setup_paths_thisDir_  = fileparts(mfilename('fullpath'));   % .../experiments
setup_paths_repoRoot_ = fileparts(setup_paths_thisDir_);    % .../lidar-2603

addpath(setup_paths_repoRoot_);
addpath(genpath(fullfile(setup_paths_repoRoot_, 'utils')));
addpath(genpath(fullfile(setup_paths_repoRoot_, 'experiments')));

clear setup_paths_thisDir_ setup_paths_repoRoot_
