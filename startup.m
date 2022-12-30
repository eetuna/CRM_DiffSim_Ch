%
% This script adds source directories to path and sets up the workspace.
%

% Add data files to path


addpath(genpath('./3D_dynamic_response_data'));

% Add matlab source to path
addpath(genpath('./matlab'));

addpath(genpath('./Numerical_quasi_static'));


% Remove the magic (clownyness stays)
clear;