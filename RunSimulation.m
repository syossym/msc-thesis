%
% This file runs the main procedures used throughout the
% calculations in the thesis
%
% Tested: Matlab 7.8.0, 7.9.0, 7.9.0
% Created by: Yossi Michaeli, October 2011
% Edited by: -
%

%% Init calculation
close all; clear all; clc;
warning off;

% Set project directories
global project_path;
project_path = 'C:\Users\Yossi Michaeli\Documents\Code\SVN\MSc_Thesis_Repo';    % PUT HERE THE ROOT DIRECTORY PATH
cd(project_path);
run('.\Common\AddPath.m');

% Create the physical constants structure
Constants;

%% Run the calculation

% Define the simulation input parameters
Params.T = 2;                                               % ambient temperature [K]
Params.Orientation = 0;                                     % the angle between the k vector components - 
                                                            % lies between 0 (100) and pi (110)
Params.ValenceCalculationMethod = 'PropagationMatrix';      % the k.p computational method type
Params.num_valence_subbands = 3;                            % number of valence subbands to consider
Params.num_cond_subbands = 1;                               % number of conduction subbands to consider
Params.gamma_vec = [1e11];                                  % phenomenological broadening factor [s^-1]

if (Params.T == 77)                                         % the k-axis quantization for the k.p calculation
    Params.k_indices = 1:100;                           
elseif (Params.T == 2)
    Params.k_indices = 1:150;
end

Params.k_scale = 5;                                         % k-axis scaling factor for the 
Params.P_vec = 1e6;                                         % hole concentration in the active region [cm^-2]

if (Params.T == 77)                                         % the spectrum calculation energy axis [J]
    Params.E_exc = 1.45*Consts.e_0 : 5e-5*Consts.e_0 : 1.5*Consts.e_0;
elseif (Params.T == 2)
    Params.E_exc = 1.5*Consts.e_0 : 5e-5*Consts.e_0 : 1.6*Consts.e_0;%1.5*Consts.e_0 : 5e-5*Consts.e_0 : 1.6*Consts.e_0;
end

Params.delta_vec = [0.96:0.002:1];                          % the MC detuning parameters vector
Params.num_mins = 4;                                        % number of minima to search in the reflection spectraq
Params.delta = 0.96; %0.95574;                              % the MC detuning for the bare MC reflection calculation
Params.N_DEG_vec = [1e11,1.5e11,2e11,2.5e11];               % electron concentration in the active region [cm^-2];
potential_vec = [0.84:0.000001:0.85];                       % the initial potential for the S-P calculation [eV]                           

% Read the simulated structure file
Structure = ReadSimulatedStructure(Params);

close all;
for (nn=1:length(Params.N_DEG_vec))                         % iterate over all electron concentrations
    Params.N_DEG = Params.N_DEG_vec(nn);
    
    if (Params.T == 2)                                      % the initial potential for the S-P calculation [eV]                                   
        Params.potential = 7.6e-13*Params.N_DEG + 0.77;                 
    elseif (Params.T == 77)
        Params.potential = 2.433e-007*Params.N_DEG^0.5106 + 0.7089;    
    end

    % Calculate the potential well profile using the S-P method
    [Bands{nn}, Structure] = CalculateLevels(Structure, Params);
    
    for (gg=1:length(Params.gamma_vec))                     % iterate over the broadening parameter vector
        
        Params.gamma = Params.gamma_vec(gg);
        Params.P = Params.P_vec(1);
       
        % Save the parameters to the Results/Temp directory
        SaveTempData(Params, Bands{nn}, 'Bands');
        SaveTempData(Params, Structure, 'Structure');
       
        % Calculate the electron properties of the bare quantum well
        close all; drawnow;
        [Bands{nn}, QWParams{nn,gg}] = CalculateLevelsDispersion(Structure, Bands{nn}, Params);
        
        % Save the parameters to the Results/Temp directory
        SaveTempData(Params, Bands{nn}, 'Bands');
        SaveTempData(Params, QWParams{nn,gg}, 'QWParams');
        
        % Calculate the optical properties of the bare qunatum well 
        close all; drawnow;
        QWParams{nn,gg} = CalculateOpticalParameters(Structure, Bands{nn}, QWParams{nn,gg}, Params);
        
        % Save the parameters to the Results/Temp directory
        SaveTempData(Params, QWParams{nn,gg}, 'QWParams');
        
        % Calculate the optical properties of the entire MC structure
        close all; drawnow;
        MCParams{nn,gg} = CalculateFullStructureReflection(Structure, Bands{nn}, QWParams{nn,gg}, Params);
        
        % Save the parameters to the Results/Temp directory
        SaveTempData(Params, MCParams{nn,gg}, 'MCParams');
    end
end
close all; drawnow;