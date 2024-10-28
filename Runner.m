%
% This is the main simulation script through which various
% calculation scenarious are being tested.
%
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, May 2009
% Edited by: -
%

%% Init calculation
close all; clear all; clc;
warning off;

% Set project directories
global project_path;
project_path = 'C:\Users\Yossi Michaeli\Documents\Code\SVN\MSc_Thesis_Repo';
%project_path = 'z:\Yossi Michaeli\Documents\Thesis\Code';
cd(project_path);
run('.\Common\AddPath.m');

% Create the physical constants structure
Constants;

%% Run the simulations
sim_num = 10;    % marks the number of the simulation to run.

switch(sim_num)
    case 1,    % Bulk SC band structure - k.p --------------------------
        
        disp('-- Sim. 1: Bulk SC band structure - k.p');
        
        material_name = 'GaAs'
        grid_size = 501
        
        Results_Sim_1 = Bulk_kp(material_name, grid_size);
        
    case 2,    % Single QW band structure - k.p ------------------------
        
        disp('-- Sim. 2: Single QW band structure - k.p');
        
        well_width = 180   % [A]
        Results_Sim_2 = SingleQW_kp(well_width, 0.1);
        
    case 3,    % Single QW absorption and gain -------------------------
        
        disp('-- Sim. 3: Single QW absorption and gain');
        
        if (exist('Results_Sim_2'))
            Results_Sim_3 = SingleQW_AbsorptionGain(Results_Sim_2);
        else
            error('It seems the single QW k.p simulation results are missing');
        end
        
    case 4,    % Single/multiple QW band structure - transfer matrix ---
        
        disp('-- Sim. 4: Single/multiple QW band structure - transfer matrix');
        
        Structure = { 'GaAlAs' , 150, 0.3 ;
            'GaAs' , 100, 0 ;
            'GaAlAs' , 150, 0.3 ;}
        T = 300;
        Params.T = T;
        Params.Orientation = pi;   % degree between the k vector components - between 0 (100) and pi (110) 
        clear Results_Sim_4;
        Results_Sim_4 = GeneralStructure_TransferMatrix(Structure,Params);

    case 5,    % DOS, absorption and gain calculation ------------------
        
        disp('-- Sim. 5: DOS, absorption and gain calculation');
        
        Results_Sim_5_Dos = CalculateDos(Results_Sim_4);
        
        GainParams.N_v = 5e1;               % valence band population (1e10 cm^-2)
        GainParams.N_c = 1e1;               % conduction band population (1e10 cm^-2)
        GainParams.Gamma = 0.01*Consts.e_0; % Lorentzian linewidth
        GainParams.T = 20;                 % temperature (K)
        
        abs_tot_TE = 0; abs_tot_TM = 0;
        for(vb_index = 1:5)
            for (cb_index = 1:2)
                
                Results_Sim_5_Gain = CalculateGain(vb_index,cb_index,2,catstruct(Results_Sim_4, Results_Sim_5_Dos, GainParams));
                abs_tot_TE = abs_tot_TE + Results_Sim_5_Gain.alpha_E_TE;
                abs_tot_TM = abs_tot_TM + Results_Sim_5_Gain.alpha_E_TM;
                
            end
        end
        
        figure;
        plot(Results_Sim_5_Gain.E_pump, abs_tot_TE, 'b', Results_Sim_5_Gain.E_pump, abs_tot_TM, 'g');
        
        %         Results_Sim_5 = DosAbsorptionGain(Results_Sim_4, GainParams);
        
    case 6,    % Optical proccesses calcultion -------------------------
        
        disp(' -- Structure simulation ------------------------------------------');
        disp(' ');
        t_start = tic;
        
        % Structure defintion
        L_z = 180e-10;
        x = 0.102;
        T = 77;
        Structure = { 'GaAlAs' , 80, x ;
                      'GaAs', L_z*1e10, 0;
                      'GaAlAs' , 80, x }
        
        % Run the simulation
        R = GeneralStructure_TransferMatrix(Structure, T);
        %close all;
        h_Band_Diagrams = PlotBandDiagram(R.k_t_mat, R.E_k_c, R.E_k_v, x, T, L_z, R.Materials{2}, 0);
        
        disp(' ');
        disp(['This stage took ', num2str(toc(t_start)/60), ' minutes']);
        disp(' ');
        
        %SingleQW_DEG_FCT_HF;
        
    case 7,    % MC reflection calculation with and without QW ---------
        
        MCQWActiveReflectionSpectrum_GeneralPurpose;
        
    case 8,    % Self-consistent Sch.-Poisson band structure calc. -----
        
        initaquila;                   % initialize AQUILA
        aquila_control.mode=1;        % 1D-Simulation
        aquila_control.fix_doping=1;  % handle doping as doping levels, not as space charge
        
        % Structure definition starting at surface
        BuildAquilaStructure();
        
        add_boundary(LEFT,POTENTIAL,0);    % surface, no gate
        add_boundary(RIGHT,POTENTIAL,0);
        %add_boundary(RIGHT,FIELD,0);       % transition to bulk
        %add_boundary(LEFT,FIELD,0);
        
        runstructure;
        
    case 9,    % Full structure simulation (1) -------------------------
        
        params.T = 2;
        params.ValenceCalculationMethod = 'PropagationMatrix';
        params.num_valence_subbands = 3;
        params.num_cond_subbands = 1;
        params.gamma = 1e11;              % [s^-1]
        params.k_indices = 1:50;
        params.k_scale = 5;
        params.P = 1e6;                   % [cm^-2]
        params.E_exc = 1.45*Consts.e_0 : 1e-4*Consts.e_0 : 1.6*Consts.e_0;
        params.delta_vec = [0.96:0.002:1];
        params.num_mins = 4;
        params.delta = 0.95574;
        
        % Read the structure file
        [QStruct, EStruct] = ReadSimulatedStructure(params.T);
        
        QStruct.N_DEG_vec = [2e9,4e9,6e9,8e9,2e10,4e10,6e10,8e10,9e10,1e11,2e11,4e11,6e11,8e11];    % [cm^-2];
        %QStruct.N_DEG_vec = 1e11;
        
        for (nn=1:length(QStruct.N_DEG_vec))
            % Calculate the active region band structure
            params.N_DEG = QStruct.N_DEG_vec(nn);     % [cm^-2]
            [Cond, Valence, QStruct] = CalculateLevels(QStruct, params);
            PlotQuantumStructure(QStruct, Cond, Valence, params, 0);
            
            % Calculate the active region subbands dispersion
            [Cond, Valence, QStruct, OpticalParams{nn}] = CalculateLevelsDispersion(QStruct, Cond, Valence, params);
            
            % Calculate the optical properties of the active region
            close all;
            OpticalParams{nn} = CalculateOpticalParameters(QStruct, Cond, Valence, OpticalParams{nn}, params);
            OpticalParams{nn}.QStruct = QStruct;
            OpticalParams{nn}.Valence = Valence;
            OpticalParams{nn}.Cond = Cond;
            OpticalParams{nn}.gamma = params.gamma;
            
            % Calculate reflection spectra of the structure
            params.N_DEG_index = nn;
            FullStructureParams{nn} = CalculateFullStructureReflection(EStruct, QStruct, Valence, Cond, OpticalParams, params);
            close all;
        end
        
        % Save workspace
        
        % Plot optical simulation results
        PlotOpticalSimulationParams(QStruct, OpticalParams, FullStructureParams, params);
        
    case 10,   % Full structure simulation (2) -------------------------
        
        Params.T = 2;
        Params.Orientation = 0;   % degree between the k vector components - lies between 0 (100) and pi (110) 
        Params.ValenceCalculationMethod = 'PropagationMatrix';
        Params.num_valence_subbands = 3;
        Params.num_cond_subbands = 1;
        Params.gamma_vec = [1e11];          % [s^-1]
        if (Params.T == 77)
            Params.k_indices = 1:100;
        elseif (Params.T == 2)
            Params.k_indices = 1:150;
        end
        Params.k_scale = 5;        
        Params.P_vec = 1e6;                 % [cm^-2]
        if (Params.T == 77)
            Params.E_exc = 1.45*Consts.e_0 : 5e-5*Consts.e_0 : 1.5*Consts.e_0;
        elseif (Params.T == 2)
            Params.E_exc = 1.5*Consts.e_0 : 5e-5*Consts.e_0 : 1.6*Consts.e_0;%1.5*Consts.e_0 : 5e-5*Consts.e_0 : 1.6*Consts.e_0;
        end
        Params.delta_vec = [0.96:0.002:1];
        Params.num_mins = 4;
        Params.delta = 0.96; %0.95574;
        Params.N_DEG_vec = [1e11,1.5e11,2e11,2.5e11];    % [cm^-2];
        potential_vec = [0.84:0.000001:0.85];
        
        Structure = ReadSimulatedStructure(Params);
        close all;
        for (nn=1:length(Params.N_DEG_vec))
        %for (nn=1:length(Params.P_vec))
            %Params.N_DEG = Params.N_DEG_vec;
            Params.N_DEG = Params.N_DEG_vec(nn);
            if (Params.T == 2)
                Params.potential = 7.6e-13*Params.N_DEG + 0.77;  % T = 2K
            elseif (Params.T == 77)
                Params.potential = 2.433e-007*Params.N_DEG^0.5106 + 0.7089; % T = 77K
            end
            %2e9 - 0.7225
            %3e9 - 0.7263
            %5e9 - 0.73115
            %8e9 - 0.7365
            %2e10 - 0.752
            %3e10 - 0.7617
            %5e10 - 0.78
            %6e10 - 0.788
            %7e10 - 0.7945
            %8e10 - 0.8
            %1e11 - 0.807
            
            %Params.potential = 0.7945;
            [Bands{nn}, Structure] = CalculateLevels(Structure, Params);
            
            for (gg=1:length(Params.gamma_vec))
                
                Params.gamma = Params.gamma_vec(gg);
                %Params.P = Params.P_vec(nn);
                Params.P = Params.P_vec(1);
                
                %                 for (pp=1:length(potential_vec))
                %                     potential_vec(pp)
                %                     Params.potential = potential_vec(pp);
                %                     [Bands{nn}, Structure] = CalculateLevels(Structure, Params);
                %                     test_vec = Bands{1}.Cond.Ve_profile-Bands{1}.Cond.Phi;
                %                     if (abs(test_vec(end)-test_vec(end-1))<1e-5)
                %                         break;
                %                     end
                %                 end
                       
                SaveTempData(Params, Bands{nn}, 'Bands');
                SaveTempData(Params, Structure, 'Structure');
                
                close all; drawnow;
                [Bands{nn}, QWParams{nn,gg}] = CalculateLevelsDispersion(Structure, Bands{nn}, Params);
                SaveTempData(Params, Bands{nn}, 'Bands');
                SaveTempData(Params, QWParams{nn,gg}, 'QWParams');
                
                close all; drawnow;
                QWParams{nn,gg} = CalculateOpticalParameters(Structure, Bands{nn}, QWParams{nn,gg}, Params);
                SaveTempData(Params, QWParams{nn,gg}, 'QWParams');
                
                close all; drawnow;
                MCParams{nn,gg} = CalculateFullStructureReflection(Structure, Bands{nn}, QWParams{nn,gg}, Params);
                SaveTempData(Params, MCParams{nn,gg}, 'MCParams');
            end
        end
        close all; drawnow;
        %PlotOpticalSimulationParams(Structure, Bands, QWParams, MCParams, Params);
        
    case 11,   % MC Relfection Calculation -----------------------------
        
        % Read the simulation file
        [filename, pathname, filterindex] = uigetfile({'*.mat', 'Mat Files'; '*.*',  'All Files'}, 'Select simulation results file');
        load([pathname filename]);
        
        Params.delta_vec = [0.945:0.00002:0.958];
             
        % Calculate the MC reflection spectra
        MCParams_New = CalculateFullStructureReflection(Structure, Bands, QWParams, Params);
        %SaveTempData(Params, MCParams, 'MCParams');
        close all;
        
    case 12,   % Simulation Results Plotting ---------------------------
        
        dirname = uigetdir(project_path, 'Pick a Directory');
        filenames = dir(dirname);
        
        type = 'Const_gamma'; %'Const_N'
        params_vec = []; filenames_vec = [];
        for (ii=1:length(filenames))
            switch (type)
                case 'Const_gamma',
                    if (~filenames(ii).isdir)
                        params_vec = [params_vec; str2double(strrep(substr(filenames(ii).name, 11, 7),'_','.'))];
                        filenames_vec = [filenames_vec; filenames(ii).name];
                    end
                case 'Const_N',
                    if (~filenames(ii).isdir)
                        params_vec = [params_vec; str2double(strrep(substr(filenames(ii).name, 25, 7),'_','.'))];
                        filenames_vec = [filenames_vec; filenames(ii).name];
                    end
            end
        end
        
        [Y,I] = sort(params_vec, 1);
        filenames_vec = filenames_vec(I,:);
        
        for (ii=1:length(filenames_vec(:,1)))
           switch (type)
                case 'Const_gamma',
                    load([dirname '\' filenames_vec(ii,:)]); 
                    if (ii==1)
                       Params_Full = Params; 
                       Params_Full.gamma_vec = Params.gamma;
                    end
                    Params_Full.N_DEG_vec(ii) = Params.N_DEG;
                    
                    Bands_Full{ii} = Bands;
                    QWParams_Full{ii,1} = QWParams;
                    MCParams_Full{ii,1} = MCParams;
                case 'Const_N',
                    load([dirname '\' filenames_vec(ii,:)]); 
                    if (ii==1)
                       Params_Full = Params; 
                       Params_Full.N_DEG_vec = Params.N_DEG;
                    end
                    Params_Full.gamma_vec(ii) = Params.gamma;
                    
                    Bands_Full{ii} = Bands;
                    QWParams_Full{1,ii} = QWParams;
                    MCParams_Full{1,ii} = MCParams;
            end
        end
        
    case 13,   % Calculate Reflection Spectra (1) ----------------------
        
        [filename, pathname, filterindex] = uigetfile({'*.mat', 'Mat Files'; '*.*',  'All Files'}, 'Select simulation results file');
        load([pathname filename]);
        
        Params_Full.delta_vec = [0.9:0.002:0.93];
        Params_Full.num_mins = 4;
        Params_Full.delta = 0.95574;
        clear MCParams_Full;
        for (nn=1:length(Params_Full.N_DEG_vec))
            Params_Full.N_DEG = Params_Full.N_DEG_vec(nn);           
            for (gg=1:length(Params_Full.gamma_vec))
                Params_Full.gamma = Params_Full.gamma_vec(gg);
                MCParams_Full{nn,gg} = CalculateFullStructureReflection(Structure, Bands_Full{nn}, QWParams_Full{nn,gg}, Params_Full);
            end
        end
        
    case 14,   % Calculate Reflection Spectra (2) ----------------------
         
        [filename, pathname, filterindex] = uigetfile({'*.mat', 'Mat Files'; '*.*',  'All Files'}, 'Select simulation results file');
        load([pathname filename]);
        
        %project_path = 'z:\Yossi Michaeli\Documents\Thesis\Code';
        Structure = ReadSimulatedStructure(Params);
        
        E_shift = 0;   % [eV]
        
        Params.delta_vec = [0.94:0.001:0.96];
        Params.num_mins = 5;
        Params.delta = 0.9863;
        Params.E_exc = Params.E_exc - E_shift*Consts.e_0;
        
        for (nn=1:length(Params.N_DEG_vec))
            Params.N_DEG = Params.N_DEG_vec(nn);           
            for (gg=1:length(Params.gamma_vec))
                Params.gamma = Params.gamma_vec(gg);
                MCParams_New{nn,gg} = CalculateFullStructureReflection(Structure, Bands{nn}, QWParams{nn,gg}, Params);
            end
        end
               
    case 15,   % Bare MC Reflection ------------------------------------
        
        Params.T = 2;
        Params.E_exc = 1.35*Consts.e_0 : 1e-4*Consts.e_0 : 1.65*Consts.e_0;
        delta_vec = 0.94;
        width_factor = 1.0;
        Structure = ReadSimulatedStructure(Params);
        MCParams.delta_vec = delta_vec;
        MCParams.T = Params.T;
        MCParams.E = Params.E_exc;
        
        close all;
        h_fig = figure('Name', 'Structure reflection');
        cavity = []; win = [];
        for (ii=1:length(delta_vec))
            
            Params.delta = delta_vec(ii);          
            MCParams.Ref{ii} = CalculateMCReflection(Structure, Params, width_factor);
            [maxs,mins] = peakdet(abs(MCParams.Ref{ii}.r_MC), 0.005, Params.E_exc/Consts.e_0);
            MCParams.Ref{ii}.mins = mins;
            [m,i] = max(diff(diff(mins(:,1))));
            MCParams.Ref{ii}.CavityMode = mins(i+2,:);
            MCParams.Ref{ii}.Win_1 = mins(i+1,:);
            cavity = [cavity; MCParams.Ref{ii}.CavityMode];
            win = [win; MCParams.Ref{ii}.Win_1];
            
            figure(h_fig);
            subplot(2,2,[1 2]); box on; hold on;
            h = plot(Params.E_exc/Consts.e_0, abs(MCParams.Ref{ii}.r_MC), 'b', ...
                     MCParams.Ref{ii}.CavityMode(:,1)*ones(1,length(MCParams.Ref{ii}.r_MC)),  abs(MCParams.Ref{ii}.r_MC), ':r', ...
                     MCParams.Ref{ii}.Win_1(:,1)*ones(1,length(MCParams.Ref{ii}.r_MC)),  abs(MCParams.Ref{ii}.r_MC), ':r');
            %set(h,'Color',padarray(dec2binvec(mod(ii,8)), [0,3-length(dec2binvec(mod(ii,8)))], 'post').');
            xlabel('E (eV)'); ylabel('Reflection');  
            title(['\delta=' num2str(delta_vec(ii))]);
            subplot(223); hold on; box on;
            plot(delta_vec(ii), mins(:,1), '.b', delta_vec(ii), MCParams.Ref{ii}.CavityMode(1), '.r', delta_vec(ii), MCParams.Ref{ii}.Win_1(1), '.g');
            xlabel('\delta'); ylabel('E_{min} (eV)');   
            subplot(224); hold on; box on;
            plot(delta_vec(ii), mins(:,2), '.b', delta_vec(ii), MCParams.Ref{ii}.CavityMode(2), '.r', delta_vec(ii), MCParams.Ref{ii}.Win_1(2), '.g');
            xlabel('\delta'); ylabel('Reflection minimun');
                       
            drawnow;
        end
        
        save([Structure.Name '_MCReflection.mat'], 'MCParams');
        
    case 16,   % MC reflection fitting ---------------------------------
        
        E_4_7_10_1_2K_Calc = load('Results\Calibration\4-7-10.1_2K_Cavity_Mode_Win_Calc_Yulia.mat');
        E_4_7_10_1_77K_Calc = load('Results\Calibration\4-7-10.1_2K_Cavity_Mode_Win_Calc.mat');
        E_4_7_10_1_300K_Calc = load('Results\Calibration\4-7-10.1_300K_Cavity_Mode_Win_Calc.mat');
        
        E_4_6_10_1_2K_Calc = load('Results\Calibration\4-6-10.1_2K_Cavity_Mode_Win_Calc_Yulia.mat');

        E_4_7_10_1_2K_Exp = load('Results\Calibration\4-7-10.1_2K_Cavity_Mode_1.mat');
        E_4_7_10_1_77K_Exp = load('Results\Calibration\4-7-10.1_77K_Cavity_Mode.mat');
        E_4_7_10_1_300K_Exp = load('Results\Calibration\4-7-10.1_300K_Cavity_Mode.mat');
        
        E_4_6_10_1_2K_Exp = load('Results\Calibration\4-6-10.1_2K_Cavity_Mode.mat');
        
        E_4_4_05_1_300K_Exp = load('Results\Calibration\4-4-05.1_300K_Cavity_Mode.mat');
        E_4_4_05_1_300K_Calc = load('Results\Calibration\4-4-05.1_300K_Cavity_Mode_Win_Calc.mat');
        
        L8_300K_Exp = load('Results\Calibration\L8_300K_Cavity_Mode.mat');
        L8_300K_Calc = load('Results\Calibration\L8_300K_Cavity_Mode_Win_Calc.mat');
                
%        position = 2:0.1:25;      % position on the wafer [mm]
%         Reflection{2}.Cavity.Exp = -4.74e-7*position.^3+2.198e-4*position.^2-3.238e-4*position+1.548;    
%         Reflection{2}.T = 77; 
%         Reflection{1}.Cavity.Exp = 2.462e-6*position.^3+1.924e-4*position.^2-2.196e-4*position+1.549;
%         Reflection{1}.T = 2;
%         Reflection{3}.Cavity.Exp = -4.742e-6*position.^3+2.629e-4*position.^2-5.359e-4*position+1.527;
%         Reflection{3}.T = 300;
%         Reflection{1}.Window.Exp = -1.531e-5*position.^3+8.02e-4*position.^2+0.003294*position+1.505;
%         Reflection{1}.T = 2;
%         Reflection{1}.Window.Calc = [E_4_7_10_1_2K.x_win, E_4_7_10_1_2K.y_win];
%         Reflection{1}.T = 2;

%         Reflection{1}.Cavity.Exp = [E_4_7_10_1_2K_Exp.x, E_4_7_10_1_2K_Exp.y];
%         Reflection{1}.T = 2;
%         Reflection{2}.Cavity.Exp = [E_4_7_10_1_77K_Exp.x, E_4_7_10_1_77K_Exp.y];
%         Reflection{2}.T = 77;
%         Reflection{3}.Cavity.Exp = [E_4_7_10_1_300K_Exp.x, E_4_7_10_1_300K_Exp.y];
%         Reflection{3}.T = 300;

%         Reflection{1}.Cavity.Calc = [E_4_7_10_1_2K_Calc.x_cavity, E_4_7_10_1_2K_Calc.y_cavity];
%         Reflection{1}.T = 2;
%         Reflection{2}.Cavity.Calc = [E_4_7_10_1_77K_Calc.x_cavity, E_4_7_10_1_77K_Calc.y_cavity];
%         Reflection{2}.T = 77;
%         Reflection{3}.Cavity.Calc = [E_4_7_10_1_300K_Calc.x_cavity, E_4_7_10_1_300K_Calc.y_cavity];
%         Reflection{3}.T = 300;
        

        Reflection{1}.Cavity.Exp = [E_4_4_05_1_300K_Exp.x, E_4_4_05_1_300K_Exp.y];
        Reflection{1}.T = 300;

        Reflection{1}.Cavity.Calc = [E_4_4_05_1_300K_Calc.x_cavity, E_4_4_05_1_300K_Calc.y_cavity];
        Reflection{1}.T = 300;
        
        E_exc = 1.4*Consts.e_0 : 1e-3*Consts.e_0 : 1.7*Consts.e_0;
        delta_vec = 0.82:0.001:1; 
        Params.E_exc = E_exc;
        [filename, pathname, filterindex] = uigetfile({'*.str', 'Structure file'; '*.*',  'All Files'}, 'Select the structure file');
        Params.filename = filename;
        Params.pathname = pathname;
        
        Mapping = FitMCReflectionToExperiment(Params, delta_vec, Reflection, 'Interp'); 
        
        figure(1);
        colors = {'b', 'r', 'g'};
        for (ii=1:length(Mapping.Cavity))
           subplot(211); box on; hold on; 
           plot(Mapping.Cavity{ii}(:,1), Mapping.Cavity{ii}(:,2), colors{ii});
           ylabel('\delta=L(r)/L(r=0)');
           subplot(212); box on; hold on; 
           plot(Mapping.Cavity{ii}(:,1), 100*Mapping.Cavity{ii}(:,2), colors{ii});
           ylabel('Cap layer width (A)'); xlabel('Distance from wafer center (mm)');
        end
        
        figure(2);
        for (ii=1:length(Reflection))
           box on; hold on; 
           plot(Reflection{ii}.Cavity.Exp(:,1), Reflection{ii}.Cavity.Exp(:,2), colors{ii});
           ylabel('E_c (eV)'); xlabel('Distance from wafer center (mm)');
        end
        
    case 17,   % Subband calculation for various well widths -----------
        
        Params.T = 2; x = 0.102; Params.Orientation = 0;
        temp = load([project_path '\4_7_10_1_Calib.mat']);
        delta = temp.g(:,2);
        width_vec = 180; %200.*delta; %[180, 200, 250]; %180*[0.86:0.01:0.96];
        E_x_vec = load([project_path '\s1_Exciton_Binding_Energy_x_0_1.mat']);
        n_e = 1e11;    % [cm^-2]
        k_F = sqrt(2*pi*n_e*100^2);
        
        for (ii=1:length(width_vec))
            width_vec(ii)
            close all;
            Structure = { 'GaAlAs' , 100, x ;
                          'GaAs' , width_vec(ii), 0 ;
                          'GaAlAs' , 100, x }
            Results_Sim_4{ii} = GeneralStructure_TransferMatrix(Structure,Params);
            E_x(ii) = 1e-3*interp1(E_x_vec.g(:,1), E_x_vec.g(:,2), width_vec(ii), 'pchip');
        end
        
        for (ii=1:length(Results_Sim_4))
           for (cb=1:2)
              for (vb=1:3) 
                E_c_0 = Results_Sim_4{ii}.E_k_c{cb}(1)/Consts.e_0;
                E_v_0 = Results_Sim_4{ii}.E_k_v{vb}(1)/Consts.e_0;
                E_c_k_F = interp1(Results_Sim_4{ii}.k_t_mat, Results_Sim_4{ii}.E_k_c{cb}, k_F, 'pchip')/Consts.e_0;
                E_v_k_F = interp1(Results_Sim_4{ii}.k_t_mat, Results_Sim_4{ii}.E_k_v{vb}, k_F, 'pchip')/Consts.e_0;

                diff_0{cb,vb}(ii) = E_c_0 + E_v_0;
                diff_k_F{cb,vb}(ii) = E_c_k_F + E_v_k_F;
              end
           end
           diff_1s(ii) = Results_Sim_4{ii}.E_k_c{1}(1)/Consts.e_0 + Results_Sim_4{ii}.E_k_v{1}(1)/Consts.e_0 + E_x(ii); 
        end
        
        close all;
        for (cb=1:2)
            for (vb=1:3)
                figure(1); box on; hold on;
                plot(width_vec, diff_0{cb,vb}, '.', ...
                     width_vec, diff_k_F{cb,vb}, 'or');
                xlabel('Well width (A)'); ylabel('\DeltaE_{i,j} (eV)');
                text(width_vec(1)-5, diff_0{cb,vb}(1), ['(' num2str(cb) ',' num2str(vb) ')']);
            end
        end
        
        for (ii=1:length(Results_Sim_4))
            PlotOpticalMatrixElements(Results_Sim_4{ii}.k_t_mat, Results_Sim_4{ii}.rtsTE, Results_Sim_4{ii}.rtsTM, x, T, width_vec(ii), 'GaAs', 0)
        end
        
    case 18,   % Subband calculation for various angles ----------------
        
        a = 200;
        Structure = { 'GaAlAs' , a+50, 0.1 ;
                      'GaAs' , a, 0 ;
                      'GaAlAs' , a+50, 0.1 ;}
                  
        cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];
                  
        T = [2];
        angles = [0,pi];
        
        for (ii=1:length(angles))
            for (jj=1:length(T))
                Params.T = T(jj);
                Params.Orientation = angles(ii);
                Res{ii,jj} = GeneralStructure_TransferMatrix(Structure,Params);
                close all; 
                Dos{ii,jj} = CalculateDos(Res{ii,jj});
                close all;
            end
        end
        
        figure(1); 
        for (ii=1:length(Res))
            for (jj=1:length(T))
                subplot(length(T),2,jj); box on; hold on;
                for (kk=1:length(Res{ii}.E_k_v))
                    plot(Res{ii}.k_t_mat(1:length(Res{ii}.E_k_v{kk}))/100, -1e3*Res{ii}.E_k_v{kk}./Consts.e_0, 'Color', cl(ii));
                end
                xlabel('k_{||} (cm^{-1})'); ylabel('E (meV)');
                subplot(length(T),2,jj+1); box on; hold on;
                temp = Dos{ii,jj}.Dos_e_qw;
                temp(temp==0) = [];
                plot(sum(Dos{ii,jj}.Dos_h_E)./temp(1)/Consts.e_0, -Dos{ii,jj}.Dos_h_E_grid*1e3, 'Color', cl(ii));
            end
        end
        %xlabel('k_{||} (cm^{-1})'); ylabel('E (meV)');
        
        figure(2);
        for (ii=1:length(Res))
            for (jj=1:length(T))
                subplot(length(T),2,jj); box on; hold on;
                for (kk=1:length(Res{ii}.E_k_c))
                    plot(Res{ii}.k_t_mat(1:length(Res{ii}.E_k_c{kk}))/100, 1e3*Res{ii}.E_k_c{kk}./Consts.e_0, 'Color', cl(ii));
                end
                xlabel('k_{||} (cm^{-1})'); ylabel('E (meV)');
                subplot(length(T),2,jj+1); box on; hold on;
                temp = Dos{ii,jj}.Dos_e_qw;
                temp(temp==0) = [];
                plot(Dos{ii,jj}.Dos_e_qw./temp(1), Dos{ii,jj}.Dos_e_E_grid*1e3+min(Res{ii}.E_k_c{1}./Consts.e_0)*1e3, 'Color', cl(ii));
                xlabel('\rho_e/\rho_e_1');
            end
        end
        
    otherwise
        disp('Undefined simulation number');
end