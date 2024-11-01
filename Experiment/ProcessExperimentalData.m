
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

%% Read the experiment al results files

global father_dirname;

father_dirname = uigetdir(project_path, 'Pick Experimental Data Directory');
filenames = dir(father_dirname);

structure_name = strsplit('\', father_dirname);
Results.structure_name = structure_name(end);

Results = ReadDirectoryContents(father_dirname, Results, 0);

%% Proccess the data

% 300K
try
    if (IsStructField(Results.T300K, 'Ref'))
        fields_ref = fieldnames(Results.T300K.Ref);
    end
    if (IsStructField(Results.T300K, 'PL'))
        fields_pl = fieldnames(Results.T300K.PL);
    end
    for(ii=1:length(fields_ref))
        if (IsStructField(Results.T300K, 'Ref'))
            Results.T300K.Ref.(fields_ref{ii}).Data = Results.T300K.Ref.(fields_ref{ii}).Data./repmat(FindLampReflection(Results.LampRef, Results.T300K.Ref.(fields_ref{ii}).Spectral_Offset).', length(Results.T300K.Ref.(fields_ref{ii}).Data(:,1)), 1);
            Results.T300K.Ref.(fields_ref{ii}) = ExtractParameters(Results.T300K.Ref.(fields_ref{ii}),Results.T300K.T, 'Ref');
        end
    end
    for(ii=1:length(fields_pl))
        if (IsStructField(Results.T300K, 'PL'))
            Results.T300K.PL.(fields_pl{ii}) = ExtractParameters(Results.T300K.PL.(fields_pl{ii}),Results.T300K.T, 'PL');
        end
    end
catch exc1
    disp(exc1.message);
end

% 77K
try
    if (IsStructField(Results.T77K, 'Ref'))
        fields_ref = fieldnames(Results.T77K.Ref);
    end
    if (IsStructField(Results.T77K, 'PL'))
        fields_pl = fieldnames(Results.T77K.PL);
    end
    for(ii=1:length(fields_ref))
        if (IsStructField(Results.T77K, 'Ref'))
            Results.T77K.Ref.(fields_ref{ii}).Data = Results.T77K.Ref.(fields_ref{ii}).Data./repmat(FindLampReflection(Results.LampRef, Results.T77K.Ref.(fields_ref{ii}).Spectral_Offset).', length(Results.T77K.Ref.(fields_ref{ii}).Data(:,1)), 1);
            Results.T77K.Ref.(fields_ref{ii}) = ExtractParameters(Results.T77K.Ref.(fields_ref{ii}),Results.T77K.T, 'Ref');
        end
    end 
    for(ii=1:length(fields_pl))
        if (IsStructField(Results.T77K, 'PL'))
            Results.T77K.PL.(fields_pl{ii}) = ExtractParameters(Results.T77K.PL.(fields_pl{ii}),Results.T77K.T, 'PL');
        end
    end
catch exc2
    disp(exc2.message);
end

% 2K
try
    if (IsStructField(Results.T2K, 'Ref'))
        fields_ref = fieldnames(Results.T2K.Ref);
    end
    if (IsStructField(Results.T2K, 'PL'))
        fields_pl = fieldnames(Results.T2K.PL);
    end
    for (ii=1:length(fields_ref))
        if (IsStructField(Results.T2K, 'Ref'))
            Results.T2K.Ref.(fields_ref{ii}).Data = Results.T2K.Ref.(fields_ref{ii}).Data;%./repmat(FindLampReflection(Results.LampRef, Results.T2K.Ref.(fields{ii}).Spectral_Offset).', length(Results.T2K.Ref.(fields{ii}).Data(:,1)), 1);
            Results.T2K.Ref.(fields_ref{ii}) = ExtractParameters(Results.T2K.Ref.(fields_ref{ii}),Results.T2K.T, 'Ref');
        end
    end
    for (ii=1:length(fields_pl))
        if (IsStructField(Results.T2K, 'PL'))
            Results.T2K.PL.(fields_pl{ii}) = ExtractParameters(Results.T2K.PL.(fields_pl{ii}),Results.T2K.T, 'PL');
        end
    end
catch exc2
    disp(exc2.message);
end

%% Plot the experimental data

PlotExperimentalResults(Results);

%SimulationExperimentCompare(Results, MCParams);
