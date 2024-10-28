function [QStruct, EStruct] = ReadSimulatedStructure(T)
%
% This function reads the simulated structure *.str file and
% creates the simulation structs.
%
%   Input:  T - ambience temperature [K].
%
%   Output: 'QStruct' - the quantum structure parameters.
%           'EStruct' - the entire structure parameters.
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

%% Read the file
[filename, pathname, filterindex] = uigetfile({'*.str', 'Structure file'; '*.*',  'All Files'}, 'Select the structure file');
imported_file = importdata([pathname filename], '\t');
disp(['Imported structure file: ' filename]);

for (mm=1:length(imported_file.data(:,1)))
    EStruct{mm}.Name = imported_file.textdata{6+mm,1};
    EStruct{mm}.x = imported_file.data(mm,1);
    EStruct{mm}.L = imported_file.data(mm,2);           % [A]
    EStruct{mm}.IsDoped = imported_file.data(mm,3);
    EStruct{mm}.IsQuantum = imported_file.data(mm,4);
    EStruct{mm}.IsActive = imported_file.data(mm,5);
end

%% Build output structures

material_params.E = 1.525;   % [eV]
material_params.T = T;       % [K]
ql_index = 1; al_index = 1;
for (mm=1:length(EStruct))
    material_params.x = EStruct{mm}.x;
    material = GetMaterial(EStruct{mm}.Name,material_params);
    EStruct{mm} = catstruct(EStruct{mm}, material);
    if (EStruct{mm}.IsQuantum)
        QLayers{ql_index} = EStruct{mm};
        ql_index = ql_index + 1;
        if (EStruct{mm}.IsActive)
            ALayers{al_index} = EStruct{mm};
            al_index = al_index + 1;
        end
    end
end
QStruct.Layers = QLayers;
QStruct.ActiveLayers = ALayers;
QStruct = AddQStructProfiles(QStruct);

function QStruct = AddQStructProfiles(QStruct)

global Consts;

%% Create the profile files

% Create a 's.r' file for the structure
delete('.\Bin\Input\s.r');
delete('.\Bin\*.r');
layer_num = length(QStruct.Layers);
s_file = zeros(layer_num,3);
for (ii=1:layer_num)
    s_file(ii,:) = [QStruct.Layers{ii}.L, QStruct.Layers{ii}.x, 0];
end
save '.\Bin\Input\s.r' 's_file' -ASCII;
copyfile('.\Bin\Input\s.r', '.\Bin\s.r');

% Generate the structure profiles file 'x.r' and the doping profile file
% 'd.r'
RunExternalCode('.\Bin\efsx.exe');
copyfile('.\Bin\x.r', '.\Bin\Output\x.r');
%copyfile('.\Bin\d.r', '.\Bin\Output\d.r');

% Generate the potential profiles 'v.r'
delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxv.exe', ' -p e');
copyfile('.\Bin\v.r', '.\Bin\Output\v_e.r');
delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxv.exe', ' -p h');
copyfile('.\Bin\v.r', '.\Bin\Output\v_h.r');
delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxv.exe', ' -p l');
copyfile('.\Bin\v.r', '.\Bin\Output\v_l.r');

% Generate effective mass profiles 'm.r'
delete('.\Bin\m.r');
RunExternalCode('.\Bin\efxm.exe', ' -p e');
copyfile('.\Bin\m.r', '.\Bin\Output\m_e.r');
delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxm.exe', ' -p h');
copyfile('.\Bin\m.r', '.\Bin\Output\m_h.r');
delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxm.exe', ' -p l');
copyfile('.\Bin\m.r', '.\Bin\Output\m_l.r');

%% Add profiles to the structure

xr_file = load('.\Bin\Output\x.r');
z_grid = xr_file(:,1);                     % [m]
x_profile = xr_file(:,2);
vr_e_file = load('.\Bin\Output\v_e.r');
Ve_profile = vr_e_file(:,2);               % [eV]
vr_h_file = load('.\Bin\Output\v_h.r');
Vh_profile = vr_h_file(:,2);               % [eV]
mr_e_file = load('.\Bin\Output\m_e.r');
m_e_profile = mr_e_file(:,2);              % [kg]
mr_h_file = load('.\Bin\Output\m_h.r');
m_h_profile = mr_h_file(:,2);              % [kg]
%d_file = load('.\Bin\Output\d.r');
%doping_profile = d_file(:,2);              % [1e18cm^-3]

g1_profile = zeros(length(z_grid),1);
g2_profile = zeros(length(z_grid),1);
g3_profile = zeros(length(z_grid),1);
g1_profile(1) = QStruct.Layers{1}.g1;
g2_profile(1) = QStruct.Layers{1}.g2;
g3_profile(1) = QStruct.Layers{1}.g3;

ii = 2; mat_index = 1; L = QStruct.Layers{mat_index}.L*1e-10;
while(ii <= length(z_grid) && mat_index <= layer_num)
    g1_profile(ii) = QStruct.Layers{mat_index}.g1;
    g2_profile(ii) = QStruct.Layers{mat_index}.g2;
    g3_profile(ii) = QStruct.Layers{mat_index}.g3;
    if (z_grid(ii) >= L)
        mat_index=mat_index+1;
        L = L + QStruct.Layers{mat_index}.L*1e-10;
    end
    ii=ii+1;
end

norm_const = Consts.hbar^2/(2*Consts.m_0);
gamma1_profile = g1_profile.*norm_const;
gamma2_profile = g2_profile.*norm_const;
gamma3_profile = g3_profile.*norm_const;

%% Add to the structure
L_clad = 100e-10;

QStruct.z_grid = z_grid;
QStruct.x_profile = x_profile;
QStruct.m_profile.e = m_e_profile;
QStruct.m_profile.h = m_h_profile;
QStruct.v_e_profile = Ve_profile;
%QStruct.v_h_profile = DecimatePotentialProfile(z_grid, L_clad, QStruct, Vh_profile.');
QStruct.v_h_profile = Vh_profile;
QStruct.g1_profile = g1_profile;
QStruct.g2_profile = g2_profile;
QStruct.g3_profile = g3_profile;
QStruct.gamma1_profile = gamma1_profile;
QStruct.gamma2_profile = gamma2_profile;
QStruct.gamma3_profile = gamma3_profile;
%QStruct.doping_profile = doping_profile;


