clear all; close all; clc;

warning off;

% Init project
global project_path;
project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
cd(project_path);
run('.\Common\AddPath.m');

% Create the constants
Constants;

% The structure
Structure = { 'GaAs' , 100, 0 ;
              'GaAlAs' , 50, 0.3 ;
              'GaAs' , 100, 0 ;
              'GaAlAs' , 50, 0.3 ;
              'GaAs' , 100, 0 ;
              'GaAlAs' , 50, 0.3 ;
              'GaAs' , 100, 0 }

% Load the physical constants
global Consts;

% Init the material parameters structure
meV = 1e-3*Consts.e_0;
[row_num,column_num] = size(Structure);
Materials = cell(row_num,1);
if (column_num == 2)
    for (ii=1:row_num)
        Materials{ii} = GetMaterial(Structure{ii,1});
        Materials{ii}.Width = Structure{ii,2};
    end
elseif (column_num == 3)
    for (ii=1:length(Structure(:,1)))
        mat_params.x = Structure{ii,3};
        Materials{ii} = GetMaterial(Structure{ii,1},mat_params);
        Materials{ii}.Width = Structure{ii,2};
    end
else
    error('TransferMatrix:Input_parameter_error',...
        'The input cell array Structure must have 2 or 3 columns');
end

% Create a 's.r' file for the structure
delete('.\Bin\Input\s.r');
delete('.\Bin\*.r');
s_file = zeros(row_num,3);
for (ii=1:row_num)
    s_file(ii,:) = [Materials{ii}.Width, Materials{ii}.x, 0];
end
save '.\Bin\Input\s.r' 's_file' -ASCII;
copyfile('.\Bin\Input\s.r', '.\Bin\s.r');

% Generate the structure profiles file 'x.r'
RunExternalCode('.\Bin\efsx.exe');
copyfile('.\Bin\x.r', '.\Bin\Output\x.r');

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

% Calculate the electronic levels and WF's using the shooting method
num_cb_subbands = 6;    % number of condiction subbands to calculate

delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxv.exe', ' -p e');
delete('.\Bin\m.r');
RunExternalCode('.\Bin\efxm.exe', ' -p e');
delete('.\Bin\Ee.r');
RunExternalCode(['.\Bin\efshoot.exe', ' -p e -s ', num2str(num_cb_subbands)]);     % calculating levels
copyfile('.\Bin\Ee.r', '.\Bin\Output\Ee.r');
RunExternalCode('.\Bin\efwf.exe');                      % calculating wf's
copyfile('.\Bin\wf_e*', '.\Bin\Output\');

% Init internal data structures
xr_file = load('.\Bin\Output\x.r');
z_grid = xr_file(:,1);          % [m]
x_profile = xr_file(:,2);
vr_e_file = load('.\Bin\Output\v_e.r');
v_e_profile = vr_e_file(:,2);   % [eV]
vr_h_file = load('.\Bin\Output\v_h.r');
v_h_profile = vr_h_file(:,2);   % [eV]
v_h_profile = v_h_profile;
mr_e_file = load('.\Bin\Output\m_e.r');
m_e_profile = mr_e_file(:,2);   % [kg]
mr_h_file = load('.\Bin\Output\m_h.r');
m_h_profile = mr_h_file(:,2);   % [kg]

g1_profile = zeros(length(z_grid),1);
g2_profile = zeros(length(z_grid),1);
g3_profile = zeros(length(z_grid),1);
g1_profile(1) = Materials{1}.g1;
g2_profile(1) = Materials{1}.g2;
g3_profile(1) = Materials{1}.g3;

ii = 2; mat_index = 1;
while(ii <= length(z_grid) && mat_index <= row_num)
    g1_profile(ii) = Materials{mat_index}.g1;
    g2_profile(ii) = Materials{mat_index}.g2;
    g3_profile(ii) = Materials{mat_index}.g3;
    if (x_profile(ii)~=x_profile(ii-1))
        mat_index=mat_index+1;
    end
    ii=ii+1;
end

norm_const = Consts.hbar^2/(2*Consts.m_0);
gamma1_profile = g1_profile.*norm_const;
gamma2_profile = g2_profile.*norm_const;
gamma3_profile = g3_profile.*norm_const;

Profile.z_grid = z_grid;
Profile.x_profile = x_profile;
Profile.v_e_profile = v_e_profile;
Profile.v_h_profile = v_h_profile;
Profile.m_e_profile = m_e_profile;
Profile.m_h_profile = m_h_profile;
Profile.g1_profile = g1_profile;
Profile.g2_profile = g2_profile;
Profile.g3_profile = g3_profile;
Profile.gamma1_profile = gamma1_profile;
Profile.gamma2_profile = gamma2_profile;
Profile.gamma3_profile = gamma3_profile;

% Simulation parameters
Del = z_grid(2) - z_grid(1);
H_m(1) = (1./m_e_profile(1)).*(-Consts.hbar^2/Del^2);
H_p(length(z_grid)) = (1./m_e_profile(length(z_grid))).*(-Consts.hbar^2/Del^2);
for (ii=1:length(z_grid))
    if (ii~=1)
        H_m(ii) = (1./(m_e_profile(ii)+m_e_profile(ii-1))).*(-Consts.hbar^2/Del^2);
    end
    if (ii~=length(z_grid))
        H_p(ii) = (1./(m_e_profile(ii)+m_e_profile(ii+1))).*(-Consts.hbar^2/Del^2);
    end
    H(ii) = v_e_profile(ii) - H_m(ii) - H_p(ii);
end

E_grid = [0:1*meV:2*max(v_e_profile)];
N = length(z_grid);

for (ee=1:length(E_grid))
    E = E_grid(ee);
    k_L = (pi/Del).*acos((E-H(1))/(2*H_p(1)))
    k_R = (pi/Del).*acos((E-H(N))/(2*H_m(N)))
    M = diag([1, H(2:N-1)-E, 1]) +...
        diag([H_m(2:N-1), -exp(1i*k_R*Del)],-1) +...
        diag([-exp(1i*k_L*Del), H_p(2:N-1)],1);
    S = [(1-exp(2*1i*k_L*Del)), zeros(1,N-1)].';
    F = linsolve(M,S);
    t(ee) = exp(-1i*k_R*Del)*F(end);
    plot(F); title(['E = ', num2str(E/Consts.e_0),' eV']);
    drawnow;
end

