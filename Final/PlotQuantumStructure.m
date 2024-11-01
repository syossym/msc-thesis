function h_structure_profile = PlotQuantumStructure(QStruct, Cond, Valence, Bands, params, save)
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
figure(h_structure_profile); 
ones_vec_1 = ones(1, length(QStruct.z_grid));
ones_vec_2 = ones(1, length(Bands.z_grid));

subplot(2,2,[1,2]); box on; hold on;
plot(Bands.z_grid, Bands.Profile.v_e_profile./Consts.e_0, 'k', 'LineWidth', 1);
%plot(QStruct.z_grid, ones_vec_1.*Cond.E_f, ':g', 'LineWidth', 1);
for (cc=1:params.num_cond_subbands)
    plot(Bands.z_grid, ones_vec_2.*(Cond.E_0(cc)-offset_energy), [cl(cc) ':'], Bands.z_grid, ones_vec_2.*(Cond.E_0(cc)-offset_energy)+(-Cond.Wf_0_Q(cc,:)./max(-Cond.Wf_0_Q(cc,:)))./20, cl(cc), 'LineWidth', 0.2);
end

%offset_energy = max(Bands.Profile.v_h_profile./Consts.e_0);

subplot(223); box on; hold on;
plot(Bands.z_grid, -Bands.Profile.v_h_profile./Consts.e_0, 'k', 'LineWidth', 1);
for (hh=1:params.num_valence_subbands)
    plot(Bands.z_grid, -ones_vec_2.*(Valence.E_k(hh,1)/Consts.e_0-offset_energy), [cl(hh) ':'], Bands.z_grid, -ones_vec_2.*(Valence.E_k(hh,1)/Consts.e_0-offset_energy)+(Valence.Wf_0_Q{hh}(:,2).'./max(max(Valence.Wf_0_Q{hh}(:,2:3).')))./80, cl(hh), 'LineWidth', 0.2);
end

subplot(224); box on; hold on;
plot(Bands.z_grid, -Bands.Profile.v_h_profile./Consts.e_0, 'k', 'LineWidth', 1);
for (hh=1:params.num_valence_subbands)
    plot(Bands.z_grid, -ones_vec_2.*(Valence.E_k(hh,1)/Consts.e_0-offset_energy), [cl(hh) ':'], Bands.z_grid, -ones_vec_2.*(Valence.E_k(hh,1)/Consts.e_0-offset_energy)+(Valence.Wf_0_Q{hh}(:,3).'./max(max(Valence.Wf_0_Q{hh}(:,2:3).')))./80, cl(hh), 'LineWidth', 0.2);
end

xlabel('z (m)'); ylabel('E (eV)');
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