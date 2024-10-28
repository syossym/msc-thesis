
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

%% Read files

dirname = uigetdir(project_path, 'Pick Experimental Data Directory');
filenames = dir(dirname);
temp = strsplit('_', dirname);
Windows.dl = str2double(temp(3));
Windows.Ref = []; Windows.PL.L1 = []; Windows.PL.L2 = [];

for (ii=1:length(filenames))
    temp1 = strsplit('_', filenames(ii).name);
    temp2 = strsplit('.', filenames(ii).name);
    if(strcmp(temp1(1), 'pl') && strcmp(temp2(end), 'dat'))
        data_pl{ii} = load([dirname '/' filenames(ii).name]);
        temp3 = strsplit('_', temp2{3});
        eval(['Windows.PL.' temp3{2} ' = [Windows.PL.' temp3{2} '; data_pl{ii}(70:910,:)];']);
    end
    if(strcmp(temp1(1), 'ref') && strcmp(temp2(end), 'dat'))
        data_ref{ii} = load([dirname '/' filenames(ii).name]);
        Windows.Ref = [Windows.Ref; data_ref{ii}(70:910,:)];
    end
end

[Windows.Ref(:,1), i] = sort(Windows.Ref(:,1));
Windows.Ref(:,2) = Windows.Ref(i,2);
[Windows.PL.L1(:,1), i] = sort(Windows.PL.L1(:,1));
Windows.PL.L1(:,2) = Windows.PL.L1(i,2);
[Windows.PL.L2(:,1), i] = sort(Windows.PL.L2(:,1));
Windows.PL.L2(:,2) = Windows.PL.L2(i,2);

figure(1);
[AX,H1,H2] = plotyy(Windows.PL.L2(:,1), Windows.PL.L2(:,2), Windows.Ref(:,1), Windows.Ref(:,2), 'plot');
hold on; plot(AX(1), Windows.PL.L1(:,1), Windows.PL.L1(:,2), 'r');
set(AX(1), 'XLim', [min(Windows.PL.L2(:,1)), max(Windows.PL.L2(:,1))], 'YTick', []);
set(AX(2), 'XLim', [min(Windows.Ref(:,1)), max(Windows.Ref(:,1))], 'YTick', []);
set(AX(1), 'YLim', [0, max(Windows.PL.L2(:,2))]);
set(AX(2), 'YLim', [0, max(Windows.Ref(:,2))]);
set(AX(1), 'YColor', 'k');
set(AX(2), 'YColor', 'k');
xlabel('E (eV)'); ylabel('PL & Reflection (a.u.)');
legend('PL (L1)', 'PL (L2)', 'Reflection');