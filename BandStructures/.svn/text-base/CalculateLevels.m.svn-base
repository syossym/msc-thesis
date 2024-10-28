function [Cond, Valence, QStruct] = CalculateLevels(QStruct, params)
%
% This function solves the Schrodinger and Poisson equations self-consistently and
% and finds the energies and wavefunction for the conduction and valence bands of the quantum
% structure. The simulation uses the QAUILA toolbox with the quantum
% structure defintion and the 2DEG concentration as an input.
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'params' - simulaion parameters.
%
%   Output: 'Cond' - structure containig the conduction band energies and
%                    wavefucntion.
%           'Valence' - structure containig the valence band energies and
%                       wavefucntion.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Run AQUILA toolbox

% Init the simulation
initaquila;
aquila_control.mode = 1;        % 1D-Simulation
%aquila_control.fix_doping = 0;  % handle doping as doping levels, not as space charge
aquila_control.T = params.T;

% Build the structure
width = 0;
for (ss=1:length(QStruct.Layers))
    if (QStruct.Layers{ss}.L > 100)
        res = 10;
    elseif (QStruct.Layers{ss}.L < 100 && QStruct.Layers{ss}.L > 50)
        res = 5;
    elseif (QStruct.Layers{ss}.L < 50 && QStruct.Layers{ss}.L > 10)
        res = 2;
    else
        res = 1;
    end

    if (ss==1 )
        boundary_potential_left = GetMaterialBandGap(QStruct.Layers{ss},params.T)-GetMaterialBandGap(QStruct.ActiveLayers{1},params.T);
    elseif (ss==length(QStruct.Layers))
        boundary_potential_right = GetMaterialBandGap(QStruct.Layers{ss},params.T)-GetMaterialBandGap(QStruct.ActiveLayers{1},params.T);
    end

    if (QStruct.Layers{ss}.IsActive)
        add_qbox([width, width+QStruct.Layers{ss}.L], 1, 4, GE+HH+LH);
        add_pbox([width, width+QStruct.Layers{ss}.L], CB+VB);
        %add_mbox(QStruct.Layers{ss}.L, 1, QStruct.Layers{ss}.x, params.N_DEG/(QStruct.Layers{ss}.L*1e-10)/100);
        add_mbox(QStruct.Layers{ss}.L, 1, QStruct.Layers{ss}.x, 0);
    else
        if (QStruct.Layers{ss}.IsDoped)
            add_mbox(QStruct.Layers{ss}.L, 1, QStruct.Layers{ss}.x, -params.N_DEG/(QStruct.Layers{ss}.L*1e-10)/100/2);
            %add_mbox(QStruct.Layers{ss}.L, 1, QStruct.Layers{ss}.x, 0);
        else
            add_mbox(QStruct.Layers{ss}.L, 1, QStruct.Layers{ss}.x, 0);
        end
    end
    width = width + QStruct.Layers{ss}.L;
end
add_pbox([0, width], CB+VB);

% Run the simulation
startpotential(0);

add_boundary(LEFT,POTENTIAL,0.8412);
add_boundary(RIGHT,POTENTIAL,0.8412);
%add_boundary(LEFT,FIELD,0);
%add_boundary(RIGHT,FIELD,0);
runstructure;
%close all;

%% Save the results

% Get wave functions
n = length(aquila_subbands.structure.xpos);
E_g_offset = -max(aquila_material.ev);
E_grid_0 = aquila_subbands.structure.xpos(1);
dz = QStruct.z_grid(2)-QStruct.z_grid(1);
z_well_grid = aquila_subbands.structure.xpos(1):dz/1e-10:aquila_subbands.structure.xpos(end);
for (cc=1:length(aquila_subbands.ge.E))
    wf_e_mat(cc, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.ge.psi((cc-1)*n+1:cc*n), z_well_grid, 'pchip'), [0, (length(QStruct.z_grid)-length(z_well_grid)-1)/2])];
end
wf_h_mat = [];
E_0_h = [];
for (hh=1:length(aquila_subbands.hh.E))
    wf_hh_mat(hh, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.hh.psi((hh-1)*n+1:hh*n), z_well_grid, 'pchip'), [0, (length(QStruct.z_grid)-length(z_well_grid)-1)/2])];
    E_0_h = [E_0_h, aquila_subbands.hh.E(hh)+E_g_offset];
    wf_h_mat = [wf_h_mat; wf_hh_mat(hh, :)];
end
for (lh=1:length(aquila_subbands.lh.E))
    wf_lh_mat(lh, :) = [0,padarray(interp1(aquila_subbands.structure.xpos, aquila_subbands.lh.psi((lh-1)*n+1:lh*n), z_well_grid, 'pchip'), [0, (length(QStruct.z_grid)-length(z_well_grid)-1)/2])];
    E_0_h = [E_0_h, aquila_subbands.lh.E(lh)+E_g_offset];
    wf_h_mat = [wf_h_mat; wf_lh_mat(lh, :)];
end
[E_0_h, I] = sort(E_0_h, 'descend');
wf_h_mat = wf_h_mat(I, :);

% Save the energies and wavefunctions
Cond.num_subbands = length(aquila_subbands.ge.E);
Cond.E_0 = sort(aquila_subbands.ge.E+E_g_offset);       % [eV]
Cond.Wf_0 = wf_e_mat;

Valence.num_subbands = length(aquila_subbands.hh.E)+length(aquila_subbands.lh.E);
Valence.num_hh_subbands = length(aquila_subbands.hh.E);
Valence.num_lh_subbands = length(aquila_subbands.lh.E);
Valence.E_0 = E_0_h;                                    % [eV]
Valence.Wf_0 = wf_h_mat;
Valence.E_hh_0 = aquila_subbands.hh.E+E_g_offset;       % [eV]
Valence.Wf_hh_0 = wf_hh_mat;
Valence.E_lh_0 = aquila_subbands.lh.E+E_g_offset;       % [eV]
Valence.Wf_lh_0 = wf_lh_mat;

QStruct.Ve_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.ec+E_g_offset, QStruct.z_grid/1e-10, 'cubic'));
QStruct.Vh_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.ev+E_g_offset, QStruct.z_grid/1e-10, 'cubic'));
QStruct.Phi = interp1(aquila_structure.xpos, phi, QStruct.z_grid/1e-10, 'pchip');
QStruct.g1_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma1, QStruct.z_grid/1e-10, 'cubic'));
QStruct.g2_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma2, QStruct.z_grid/1e-10, 'cubic'));
QStruct.g3_profile = NormalizeStepFunction(interp1(aquila_structure.xpos, aquila_material.Gamma3, QStruct.z_grid/1e-10, 'cubic'));

offset = max(QStruct.Vh_profile-QStruct.Phi);

QStruct.Phi = QStruct.Phi + offset;
QStruct.E_f = aquila_control.Efermi + E_g_offset - offset;
Valence.E_0 = Valence.E_0 - offset;
Cond.E_0 = Cond.E_0 - offset;

norm_const = Consts.hbar^2/(2*Consts.m_0);
QStruct.gamma1_profile = QStruct.g1_profile.*norm_const;
QStruct.gamma2_profile = QStruct.g2_profile.*norm_const;
QStruct.gamma3_profile = QStruct.g3_profile.*norm_const;

function F = NormalizeStepFunction(F)

diff_F = diff(abs(F));

for (ii=2:length(diff_F))
    if (diff_F(ii-1)==0 && diff_F(ii)>0)
        F(ii+1) = F(ii+2);
    elseif (diff_F(ii-1)==0 && diff_F(ii)<0)
        F(ii+1) = F(ii+2);
    end
end
