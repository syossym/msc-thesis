function [Bands, Structure] = CalculateLevels(Structure, Params)
%
% This function solves the Schrodinger and Poisson equations self-consistently and
% and finds the energies and wavefunction for the conduction and valence bands of the quantum
% structure. The simulation uses the QAUILA toolbox with the quantum
% structure defintion and the 2DEG concentration as an input.
%
%   Input: 'Structure' - the simulated structure parameters.
%          'Params' - simulaion parameters.
%
%   Output: 'Bands' - structure containig the conduction and valence band energies and
%                     wavefucntion.
%           'Structure' - the simulated structure parameters.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Run AQUILA toolbox

% Init the simulation
initaquila;
aquila_control.mode = 1;         % 1D-Simulation
%aquila_control.fix_doping = 0;  % handle doping as doping levels, not as space charge
aquila_control.T = Params.T;
z_offset = 0;

% Build the structure
width = 0;
for (ss=1:length(Structure.QuantumStruct.Layers))
    current_width = Structure.QuantumStruct.Layers{ss}.L;
    
    if (Params.N_DEG>1e11)
        if (ss==1 || ss==length(Structure.QuantumStruct.Layers))
            current_width = current_width*100;
        end
        if (ss==1)
           z_offset = current_width - Structure.QuantumStruct.Layers{ss}.L;
        end
    end
    
    if (current_width > 100)
        res = 10;
    elseif (current_width < 100 && current_width > 50)
        res = 5;
    elseif (current_width < 50 && current_width > 10)
        res = 2;
    else
        res = 0.1;
    end
    
    if (ss==1 )
        boundary_potential_left = GetMaterialBandGap(Structure.QuantumStruct.Layers{ss}.x,Params.T)-GetMaterialBandGap(Structure.QuantumStruct.ActiveLayers{1}.x,Params.T);
    elseif (ss==length(Structure.QuantumStruct.Layers))
        boundary_potential_right = GetMaterialBandGap(Structure.QuantumStruct.Layers{ss}.x,Params.T)-GetMaterialBandGap(Structure.QuantumStruct.ActiveLayers{1}.x,Params.T);
    end
    
    if (Structure.QuantumStruct.Layers{ss}.IsActive)
        add_qbox([width, width+current_width], 1, 4, GE+HH+LH);
        add_pbox([width, width+current_width], CB+VB);
        %add_mbox(current_width, 1, Structure.QuantumStruct.Layers{ss}.x, Params.N_DEG/(current_width*1e-10)/100);
        add_mbox(current_width, 1, Structure.QuantumStruct.Layers{ss}.x, 0);
    else
        if (Structure.QuantumStruct.Layers{ss}.IsDoped)
            add_mbox(current_width, 1, Structure.QuantumStruct.Layers{ss}.x, -Params.N_DEG/(current_width*1e-10)/100/2);
            %add_mbox(current_width, 1, Structure.QuantumStruct.Layers{ss}.x, 0);
        else
            add_mbox(current_width, 1, Structure.QuantumStruct.Layers{ss}.x, 0);
        end
    end
    width = width + current_width;
end
add_pbox([0, width], CB+VB);

% Run the simulation
%startpotential(0);

% N_2DEG = 2e9cm^-2 - 0.773385
% N_2DEG = 3e9cm^-2 - 0.7742
% N_2DEG = 4e9cm^-2 - 0.775
% N_2DEG = 5e9cm^-2 - 0.7758
% N_2DEG = 6e9cm^-2 - 0.7765
% N_2DEG = 7e9cm^-2 - 0.77735
% N_2DEG = 8e9cm^-2 - 0.77815
% N_2DEG = 9e9cm^-2 - 0.778903
% N_2DEG = 1e10cm^-2 - 0.77968
% N_2DEG = 2e10cm^-2 - 0.78735
% N_2DEG = 3e10cm^-2 - 0.79502
% N_2DEG = 4e10cm^-2 - 0.802682
% N_2DEG = 5e10cm^-2 - 0.81035
% N_2DEG = 6e10cm^-2 - 0.818
% N_2DEG = 7e10cm^-2 - 0.82567
% N_2DEG = 7.4e10cm^-2 - 0.83
% N_2DEG = 8e10cm^-2 - 0.83334
% N_2DEG = 9e10cm^-2 - 0.841
% N_2DEG = 1e11cm^-2 - 0.845
% N_2DEG = 1.9e11cm^-2 - 0.8558
% N_2DEG = 2e11cm^-2 - 0.8565
% N_2DEG = 3e11cm^-2 - 0.861
% N_2DEG = 4e11cm^-2 - 0.865
% N_2DEG = 5e11cm^-2 - 0.868

add_boundary(LEFT,POTENTIAL,Params.potential);
add_boundary(RIGHT,POTENTIAL,Params.potential);
%add_boundary(LEFT,FIELD,0);
%add_boundary(RIGHT,FIELD,0);
runstructure;
%close all;

%% Save the results

% Get wave functions
n = length(aquila_subbands.structure.xpos);
E_g_offset = -max(aquila_material.ev);
E_grid_0 = aquila_subbands.structure.xpos(1);
dz = Structure.QuantumStruct.z_grid(2)-Structure.QuantumStruct.z_grid(1);
z_well_grid = aquila_subbands.structure.xpos(1):dz/1e-10:aquila_subbands.structure.xpos(end);
for (cc=1:length(aquila_subbands.ge.E))
    wf_e_mat(cc, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.ge.psi((cc-1)*n+1:cc*n), z_well_grid, 'pchip'), [0, round((length(Structure.QuantumStruct.z_grid)-length(z_well_grid)-1)/2)])];
end
wf_h_mat = [];
E_0_h = [];
for (hh=1:length(aquila_subbands.hh.E))
    wf_hh_mat(hh, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.hh.psi((hh-1)*n+1:hh*n), z_well_grid, 'pchip'), [0, round((length(Structure.QuantumStruct.z_grid)-length(z_well_grid)-1)/2)])];
    E_0_h = [E_0_h, aquila_subbands.hh.E(hh)+E_g_offset];
    wf_h_mat = [wf_h_mat; wf_hh_mat(hh, :)];
end
for (lh=1:length(aquila_subbands.lh.E))
    wf_lh_mat(lh, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.lh.psi((lh-1)*n+1:lh*n), z_well_grid, 'pchip'), [0, round((length(Structure.QuantumStruct.z_grid)-length(z_well_grid)-1)/2)])];
    E_0_h = [E_0_h, aquila_subbands.lh.E(lh)+E_g_offset];
    wf_h_mat = [wf_h_mat; wf_lh_mat(lh, :)];
end
[E_0_h, I] = sort(E_0_h, 'descend');
wf_h_mat = wf_h_mat(I, :);

% Save the energies and wavefunctions
Bands.Cond.num_subbands = length(aquila_subbands.ge.E);
Bands.Cond.E_0 = sort(aquila_subbands.ge.E+E_g_offset);       % [eV]
Bands.Cond.Wf_0 = wf_e_mat(:,1:length(Structure.QuantumStruct.z_grid));

Bands.Valence.num_subbands = length(aquila_subbands.hh.E)+length(aquila_subbands.lh.E);
Bands.Valence.num_hh_subbands = length(aquila_subbands.hh.E);
Bands.Valence.num_lh_subbands = length(aquila_subbands.lh.E);
Bands.Valence.E_0 = E_0_h;                                    % [eV]
Bands.Valence.Wf_0 = wf_h_mat(:,1:length(Structure.QuantumStruct.z_grid));
Bands.Valence.E_hh_0 = aquila_subbands.hh.E+E_g_offset;       % [eV]
Bands.Valence.Wf_hh_0 = wf_hh_mat;
Bands.Valence.E_lh_0 = aquila_subbands.lh.E+E_g_offset;       % [eV]
Bands.Valence.Wf_lh_0 = wf_lh_mat;

Bands.Cond.Ve_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.ec+E_g_offset, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'cubic'));
Bands.Valence.Vh_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.ev+E_g_offset, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'cubic'));
Phi = interp1(aquila_structure.xpos, phi, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'pchip');
Structure.QuantumStruct.g1_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma1, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'cubic'));
Structure.QuantumStruct.g2_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma2, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'cubic'));
Structure.QuantumStruct.g3_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma3, Structure.QuantumStruct.z_grid/1e-10 + z_offset, 'cubic'));

offset = max(Bands.Valence.Vh_profile-Phi);

Bands.Cond.Phi = Phi + offset;
Bands.Cond.E_f = aquila_control.Efermi + E_g_offset - offset;
Bands.Valence.Phi = Bands.Cond.Phi;
Bands.Valence.E_f = Bands.Cond.E_f;
Bands.Valence.E_0 = Bands.Valence.E_0 - offset;
Bands.Cond.E_0 = Bands.Cond.E_0 - offset;

norm_const = Consts.hbar^2/(2*Consts.m_0);
Structure.QuantumStruct.gamma1_profile = Structure.QuantumStruct.g1_profile.*norm_const;
Structure.QuantumStruct.gamma2_profile = Structure.QuantumStruct.g2_profile.*norm_const;
Structure.QuantumStruct.gamma3_profile = Structure.QuantumStruct.g3_profile.*norm_const;

%% Perform Additional Level Calculations

% % Calculating conduction subbnads using the Harrison codes
% delete('.\Bin\v.r');
% v_e_save = [Structure.QuantumStruct.z_grid, (Bands.Cond.Ve_profile-Bands.Cond.Phi)*Consts.e_0];
% save '.\Bin\v.r' 'v_e_save' -ASCII;
% delete('.\Bin\m.r');
% RunExternalCode('.\Bin\efxm.exe', ' -p e');
% delete('.\Bin\Ee.r');
% RunExternalCode(['.\Bin\efshoot.exe', ' -p e -s ', num2str(6)]);    % calculating levels
% copyfile('.\Bin\Ee.r', '.\Bin\Output\Ee.r');
% delete('.\Bin\wf_e*');
% RunExternalCode('.\Bin\efwf.exe');             % calculating wf's
% copyfile('.\Bin\wf_e*', '.\Bin\Output\');
% copyfile('.\Bin\v.r', '.\Bin\Output\v_e.r');
% 
% % Save the results
% temp_E = load('.\Bin\Output\Ee.r')           
% Bands.Cond.E_0 = temp_E(:,2)/1e3;           % [eV]

function F = NormalizeStepFunction(F)

diff_F = diff(abs(F));

for (ii=2:length(diff_F))
    if (diff_F(ii-1)==0 && diff_F(ii)>0)
        F(ii+1) = F(ii+2);
    elseif (diff_F(ii-1)==0 && diff_F(ii)<0)
        F(ii+1) = F(ii+2);
    end
end
