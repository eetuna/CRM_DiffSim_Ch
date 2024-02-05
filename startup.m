%
% This script adds source directories to path and sets up the workspace.
%

% Add data files to path
addpath(genpath('./3D_dynamic_response_data'));
addpath(genpath('./parameters'));
addpath(genpath('./catheterdata'));
addpath(genpath('./output_data'));

% Add matlab source to path
addpath(genpath('./matlab'));
addpath(genpath('./Mexfiles'));

% Remove the magic (clownyness stays)
clear;