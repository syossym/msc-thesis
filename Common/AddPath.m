%									   
% This file adds the matlabpath the directories of the project.
%
%   Input: - 
%   
%   Output: -
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

global project_path;

addpath project_path;
addpath([project_path '\Common']);
addpath([project_path '\BandStructures']);
addpath([project_path '\Excitons']);
addpath([project_path '\OpticalProccesses']);
addpath([project_path '\Polaritons']);
addpath([project_path '\Final']);
addpath([project_path '\Utilities']);
addpath([project_path '\Experiment']);

run addutils; % adds the utilities directories to the matalabpath

% Add AQUILA toolbox directory
addpath([project_path '\Misc\Toolboxes\AQUILA']);