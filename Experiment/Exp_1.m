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
temp = strsplit('\', dirname);
dl_vec = []; Windows.Ref.t = 'temp'; Windows.PL.L1.t = 'temp'; Windows.PL.L2.t = 'temp';
for (ii=1:length(filenames))
    temp1 = strsplit('_', filenames(ii).name);
    temp2 = strsplit('.', filenames(ii).name);
    temp3 = strsplit('.', temp1{end});
    dl = [temp3{1},'_',temp3{2}];
    
    if(strcmp(temp1(1), 'ref') && strcmp(temp2(end), 'dat'))
        if(~IsStructField(Windows.Ref, ['dl' dl]))
            eval(['Windows.Ref.dl' dl ' = [];']);
            dl_vec = [dl_vec, str2double([temp3{1},'.',temp3{2}])];
        end
        data_ref{ii} = load([dirname '/' filenames(ii).name]);
        data_ref{ii}(:,2) = data_ref{ii}(:,2)./max(data_ref{ii}(:,2));
        eval(['Windows.Ref.dl' dl ' = [Windows.Ref.dl' dl '; data_ref{ii}(70:910,:)];']);
    end
    if(strcmp(temp1(1), 'pl') && strcmp(temp2(end), 'dat'))
        name_dl = ['dl' dl];
        eval(['if(~IsStructField(Windows.PL.' temp{end} ', name_dl)), eval([''Windows.PL.'' temp{end} ''.dl'' dl '' = [];'']); dl_vec = [dl_vec, str2double([temp3{1},''.'',temp3{2}])]; end']);
        data_pl{ii} = load([dirname '/' filenames(ii).name]);
        eval(['Windows.PL.' temp{end} '.dl' dl ' = [Windows.PL.' temp{end} '.dl' dl '; data_pl{ii}(70:910,:)];']);
    end
end

%% Plotting

figure(1); box on; hold on;
[ss,ss_i] = sort(dl_vec, 'descend');
fields = fieldnames(Windows.Ref);
%eval(['fields = fieldnames(Windows.PL.' temp{end} ');']);
fields = fields(2:end);
fields = fields(ss_i);
delta_y = 0;
for (ii=1:length(fields))
    eval(['[Windows.Ref.' fields{ii} '(:,1), i] = sort(Windows.Ref.' fields{ii} '(:,1));']);
    eval(['Windows.Ref.' fields{ii} '(:,2) = Windows.Ref.' fields{ii} '(i,2);']);
    eval(['plot(Windows.Ref.' fields{ii} '(:,1), Windows.Ref.' fields{ii} '(:,2) + delta_y, ''.'', ''MarkerSize'', 0.5, ''Color'', [0.5 0.5 0.5]);']);
%     eval(['[Windows.PL.' temp{end} '.' fields{ii} '(:,1), i] = sort(Windows.PL.' temp{end} '.' fields{ii} '(:,1));']);
%     eval(['Windows.PL.' temp{end} '.' fields{ii} '(:,2) = Windows.PL.' temp{end} '.' fields{ii} '(i,2);']);
%     eval(['plot(Windows.PL.' temp{end} '.' fields{ii} '(:,1), Windows.PL.' temp{end} '.' fields{ii} '(:,2) + delta_y, ''-'', ''MarkerSize'', 0.5, ''Color'', [0.5 0.5 0.5]);']);
    delta_y = delta_y + 0.5;
end
xlabel('E (eV)'); ylabel('Reflection (a.u.)'); eval(['ylabel(''PL (' temp{end} ') (a.u.)'');']); 
set(gca, 'XLim', [1.49, 1.62], 'YTick', []);