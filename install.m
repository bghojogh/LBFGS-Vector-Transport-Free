function install_toolbox

clc
clear all
close all

fprintf('Installing the toolbox...\n');

% Find required paths
this_file = mfilename('fullpath'); % get this file's path and name (without .m)
toolbox_root_path = fileparts(this_file);

% Add function folders to search path
fprintf('Adding MixEst folders to search path...\n');
addpath(toolbox_root_path)
addpath(genpath(fullfile(toolbox_root_path, 'solvers')))
addpath(genpath(fullfile(toolbox_root_path, 'manifolds')))
addpath(genpath(fullfile(toolbox_root_path, 'mixest')))
addpath(genpath(fullfile(toolbox_root_path, 'thirdparty')))
addpath(genpath(fullfile(toolbox_root_path, 'utils')))
addpath(genpath(fullfile(toolbox_root_path, 'examples')))
fprintf('[Done]\n')

