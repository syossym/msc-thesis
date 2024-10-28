%% Init

close all; clear all; clc;
warning off;

% Set project directories
global project_path;
project_path = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code';
cd(project_path);
run('.\Common\AddPath.m');

% Create the physical constants structure
Constants;

colors = ['b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm'];

%% Read the experimental results files

load('Ne_White_Lamp_Spectrum.mat');
partial = [808.27, 811.9, 812.95, 813.69; 0.86538, 1.1/5.2, 0.4/5.2, 1];
partial_E(1,:) = Consts.h_J*Consts.c./(partial(1,:)*1e-9)/Consts.e_0;

global father_dirname;

father_dirname = uigetdir(project_path, 'Pick Experimental Data Directory');
filenames = dir(father_dirname);

 %hold on; box on;


for (ii=1:length(filenames))
    %stem(Ne_peak_lambda, 7e4*ones(size(Ne_peak_lambda)), ':');
    %stem(partial, 7e4*ones(size(partial)), ':r');
    temp1 = strsplit('_', filenames(ii).name);
    temp2 = strsplit('.', filenames(ii).name);
    if(strcmp(temp1(1), 'calib') && strcmp(temp2(end), 'dat'))
        data{ii} = load([father_dirname '/' filenames(ii).name]);
        [maxs,mins] = peakdet(data{ii}(:,2),50, data{ii}(:,1));
        figure(ii);
        hold on;
        stem(partial(1,:), partial(2,:), ':r');
        plot(Consts.h_J*Consts.c./(data{ii}(:,1)*Consts.e_0)/1e-9, data{ii}(:,2)./max(data{ii}(:,2)), 'Color', colors(ii));
        %stem(maxs(:,1), 7e4*ones(size(maxs(:,1))), 'r');
        hold off;
    end
end



axis tight;
