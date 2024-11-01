function h_structure_profile = PlotQuantumStructure(QStruct, Cond, Valence, params, save)
%
% This function plots the simulated quantum structure, band edge potential,
% the calculated energy levels and suitable wave functions for the
% conduction and valence bands,
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'Cond' - calculated conduction band parameters.
%          'Valence' - calculated valence band parameters.
%          'params' - simulaion parameters.
%          'save' - flag that indicates whether the plot shoud be saved.
%
%   Output: 'h_structure_profile' - the generated figure handler.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Figures init

h_structure_profile = figure('Name','Quantum Structure Profile');
cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];

%% Plotting

offset_energy = 0;
figure(h_structure_profile); box on; hold on;
ones_vec = ones(1, length(QStruct.z_grid));
plot(QStruct.z_grid, QStruct.Ve_profile, ':k', QStruct.z_grid, QStruct.Ve_profile-QStruct.Phi-offset_energy, 'k', 'LineWidth', 1);
plot(QStruct.z_grid, QStruct.Vh_profile, ':k', QStruct.z_grid, QStruct.Vh_profile-QStruct.Phi-offset_energy, 'k', 'LineWidth', 1);
plot(QStruct.z_grid, ones_vec.*QStruct.E_f, ':g', 'LineWidth', 1);
for (cc=1:length(Cond.E_0))
    plot(QStruct.z_grid, ones_vec.*(Cond.E_0(cc)-offset_energy), [cl(cc) ':'], QStruct.z_grid, ones_vec.*(Cond.E_0(cc)-offset_energy)+3.*abs(Cond.Wf_0(cc,:)).^2, cl(cc), 'LineWidth', 0.2);
end
for (hh=1:length(Valence.E_0))
    plot(QStruct.z_grid, ones_vec.*(Valence.E_0(hh)-offset_energy), [cl(hh) ':'], QStruct.z_grid, ones_vec.*(Valence.E_0(hh)-offset_energy)+3.*abs(Valence.Wf_0(hh,:)).^2, cl(hh), 'LineWidth', 0.2);
end

xlabel('z [m]'); ylabel('E [eV]');
title(['T=' num2str(params.T) 'K, N_{2DEG}=' strrep([num2str(params.N_DEG, '%1.0e')], 'e+0', '\times10^{') '}cm^{-2}']);
drawnow;

%% Optional figure saving

if (save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_structure_profile, [dirPath '\StructureProfile'], 'fig');
        exportfig(h_structure_profile, [dirPath '\StructureProfile'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
    end
end