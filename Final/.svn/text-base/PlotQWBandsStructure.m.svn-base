function PlotQWBandsStructure(Structure, Bands, QWParams, Params)
%
% This function plots the simulated quantum structure, band edge potentiala,
% the calculated energy levels and suitable wave functions for the
% conduction and valence bands.
%
%   Input: 'Structure' - the simulated structure.
%          'Bands' - the quntum band structure.
%          'Params' - simulaion parameters.
%
%   Output: 'h_structure_profile' - the generated figure handler.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Figures init

cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];
h_Band_Edges = figure('Name', 'Band Edged');
h_Subbands = figure('Name', 'Subbands');
h_Fermi = figure('Name', 'Fermi Level');
h_matrix_elements = figure('Name', 'Optical Matrix Elements');

%% Plotting

for (ii=1:length(Bands))
    Params.N_DEG_vec(ii)
    Cond = Bands{ii}.Cond;
    Valence = Bands{ii}.Valence;
    k_t_vec = Bands{ii}.k_t_vec;
    z_grid = Structure.QuantumStruct.z_grid;
    offset_energy = 0;
    ones_vec = ones(1, length(z_grid));
    
    %     figure('Name', ['Quantum Structure Profile - N2DEG=' num2str(Params.N_DEG_vec(ii), '%1.1e') 'cm^-2']);
    %     subplot(211); box on; hold on;
    %     plot(z_grid, Cond.Ve_profile, ':k', z_grid, Cond.Ve_profile-Cond.Phi-offset_energy, 'k', 'LineWidth', 1);
    %     plot(z_grid, ones_vec.*Cond.E_f, ':g', 'LineWidth', 1);
    %     text(z_grid(1), Cond.E_f, ['E_f=' num2str(Cond.E_f) 'eV']);
    %     for (cc=1:length(Cond.Wf_0(:,1)))
    %         plot(z_grid, ones_vec.*(Cond.E_0(cc)-offset_energy), [cl(cc) ':'], z_grid, ones_vec.*(Cond.E_0(cc)-offset_energy)+3.*abs(Cond.Wf_0(cc,:)).^2, cl(cc), 'LineWidth', 0.2);
    %     end
    %     axis([min(z_grid) max(z_grid) 1.5 max(Cond.Ve_profile)+0.02]);
    %     set(gca, 'XTickLabel', '');
    %     ylabel('E [eV]');
    %     title(['T=' num2str(Params.T) 'K, N_{2DEG}=' strrep([num2str(Params.N_DEG, '%1.1e')], 'e+0', '\times10^{') '}cm^{-2}']);
    %     subplot(212); box on; hold on;
    %     plot(z_grid, Valence.Vh_profile, ':k', z_grid, Valence.Vh_profile-Valence.Phi-offset_energy, 'k', 'LineWidth', 1);
    %     for (hh=1:length(Valence.E_k(:,1)))
    %         plot(z_grid, ones_vec.*(-Valence.E_k(hh,1)/Consts.e_0-offset_energy), [cl(hh) ':'],z_grid, ones_vec.*(-Valence.E_k(hh,1)/Consts.e_0-offset_energy)+3.*abs(Valence.Wf_0(hh,:)).^2, cl(hh), 'LineWidth', 0.2);
    %     end
    %     axis([min(z_grid) max(z_grid) min(Valence.Vh_profile-Valence.Phi-offset_energy) 0.05]);
    %     xlabel('z [m]'); ylabel('E [eV]');
    %     drawnow;
    
    %k_f = interp1(Cond.E_k/Consts.e_0, k_t_vec, Cond.E_f, 'pchip');
    if (IsStructField(QWParams{ii}, 'E_fc'))
        k_f = interp1(Cond.E_k(Params.k_indices)/Consts.e_0, k_t_vec, QWParams{ii}.E_fc/Consts.e_0, 'pchip');
        %E_f_vec(ii) = Cond.E_f-Cond.E_k(1)/Consts.e_0;
        E_f_vec(ii) = QWParams{ii}.E_fc/Consts.e_0-Cond.E_k(1)/Consts.e_0;
        k_f_vec(ii) = k_f;
    else
        k_f = 0;
        E_f_vec(ii) = 0;
        k_f_vec(ii) = 0;
    end
    index = 1;
    for (ee=1:length(Cond.E_k(:,1)))
        for (hh=1:length(Valence.E_k(:,1)))
            delta_E_0 = (Cond.E_k(ee,1)+Valence.E_k(hh,1))/Consts.e_0;
            delta_E_k_f = (interp1(k_t_vec, Cond.E_k(ee,Params.k_indices)/Consts.e_0, k_f, 'pchip')+interp1(k_t_vec, Valence.E_k(hh,Params.k_indices)/Consts.e_0, k_f, 'pchip'));
            details_text(index,:) = {['\DeltaE_{' num2str(ee) 'e:' num2str(hh) 'h}(k=0)=' num2str(delta_E_0) 'eV']};
            details_text(index+1,:) = {['\DeltaE_{' num2str(ee) 'e:' num2str(hh) 'h}(k=k_f)=' num2str(delta_E_k_f) 'eV']};
            index = index + 2;
        end
    end
    
    figure('Name', ['Condunction & Valence Bands - N2DEG=' num2str(Params.N_DEG_vec(ii), '%1.1e') 'cm^-2']);
    subplot(211); hold on; box on;
    plot(k_t_vec/100, 1e3*Cond.E_k(:,Params.k_indices)/Consts.e_0);
    plot(k_t_vec/100, 1e3*Cond.E_f*ones(1, length(k_t_vec)), ':r');
    text(k_t_vec(1)/100, 1e3*Cond.E_f+0.01, ['E_f=' num2str(Cond.E_f) 'eV']);
    plot(k_f*ones(1, length(k_t_vec))/100, 1e3*Cond.E_k(:,Params.k_indices)/Consts.e_0, 'k:');
    ylabel('E [meV]'); set(gca, 'XTickLabel', '');
    axis([min(k_t_vec/100) 2e6 1520 1550]);
    title(['T=' num2str(Params.T) 'K, N_{2DEG}=' strrep([num2str(Params.N_DEG_vec(ii), '%1.1e')], 'e+0', '\times10^{') '}cm^{-2}']);
    text(min(k_t_vec)/100, 1550, details_text, 'HorizontalAlignment','left','interpreter','tex');
    for (ee=1:length(Cond.E_k(:,1)))
        text(k_t_vec(10)/100, 1e3*Cond.E_k(ee,10)/Consts.e_0, [num2str(ee) 'e'], 'FontSize', 12);
    end
    subplot(212); hold on; box on;
    plot(k_t_vec/100, -1e3*Valence.E_k/Consts.e_0);
    plot(k_f*ones(1, length(k_t_vec))/100, -1e3*Valence.E_k(1,:)/Consts.e_0, 'k:');
    text(k_f/100, -23, ['k_f=' strrep([num2str(k_f/100, '%1.2e')], 'e+0', '\times10^{') '}cm^{-1}']);
    xlabel('k_{||} [cm^{-1}]'); ylabel('E [meV]');
    axis([min(k_t_vec/100) 2e6 -25 0]);
    for (hh=1:length(Valence.E_k(:,1)))
        text(k_t_vec(10)/100, -1e3*Valence.E_k(hh,10)/Consts.e_0, [num2str(hh) 'h'], 'FontSize', 12);
    end
    drawnow;
    
    figure(h_Band_Edges);
    subplot(211); box on; hold on;
    plot(z_grid/1e-10, Cond.Ve_profile, '--r', z_grid/1e-10, Cond.Ve_profile-Cond.Phi-offset_energy, 'k', 'LineWidth', 0.2);
    axis([min(z_grid/1e-10) max(z_grid/1e-10) 1.5 max(Cond.Ve_profile)+0.02]);
    set(gca, 'XTickLabel', '');
    ylabel('E (eV)');
    subplot(212); box on; hold on;
    plot(z_grid/1e-10, Valence.Vh_profile, '--r', z_grid/1e-10, Valence.Vh_profile-Valence.Phi-offset_energy, 'k', 'LineWidth', 0.2);
    axis([min(z_grid/1e-10) max(z_grid/1e-10) min(Valence.Vh_profile-Valence.Phi-offset_energy) 0.05]);
    xlabel('z (Angstrom)');
    drawnow;
    
    figure(h_Subbands);
    subplot(211); hold on; box on;
    plot(k_t_vec/100, 1e3*Cond.E_k(:,Params.k_indices)/Consts.e_0, 'LineWidth', 0.2);
    ylabel('E (meV)'); set(gca, 'XTickLabel', '');
    axis([min(k_t_vec/100) 2e6 1520 1550]);
    subplot(212); hold on; box on;
    plot(k_t_vec/100, -1e3*Valence.E_k(:,Params.k_indices)/Consts.e_0, 'LineWidth', 0.2);
    xlabel('k_t (cm^{-1})');
    axis([min(k_t_vec/100) 2e6 -25 0]);
    drawnow;
    
    [vb_len, cb_len, k_len] = size(QWParams{ii}.OptMatrixElement.TE);
    
    figure(h_matrix_elements);
    subplot(211); box on; hold on;
    for (cb_num=1:cb_len)
        plot(k_t_vec/100, squeeze(QWParams{ii}.OptMatrixElement.TE(:,cb_num,:)), 'LineWidth', 0.2);
        for (vb_num=1:vb_len)
            [m,m_i] = max(squeeze(QWParams{ii}.OptMatrixElement.TE(vb_num,cb_num,:)));
            if (ii==1)
                text(k_t_vec(m_i)/100, m, ['e' num2str(cb_num) '-h' num2str(vb_num)]);
            end
        end
    end
    title('(a)');
    ylabel('|\mu_{k}^{TE}|^2/|\mu_{k_{0}}|^2');
    set(gca, 'XTickLabel', '');
    subplot(212); box on; hold on;
    for (cb_num=1:cb_len)
        plot(k_t_vec/100, squeeze(QWParams{ii}.OptMatrixElement.TM(:,cb_num,:)), 'LineWidth', 0.2);
        for (vb_num=1:vb_len)
            [m,m_i] = max(squeeze(QWParams{ii}.OptMatrixElement.TM(vb_num,cb_num,:)));
            if (ii==1)
                text(k_t_vec(m_i)/100, m, ['e' num2str(cb_num) '-h' num2str(vb_num)]);
            end
        end
    end
    title('(b)');
    xlabel('k_t (cm^{-1})');
    ylabel('|\mu_{k}^{TM}|^2/|\mu_{k_{0}}|^2');
    
end

figure(h_Fermi);
subplot(211);
plot(Params.N_DEG_vec, E_f_vec*1e3, 'b', Params.N_DEG_vec, 3.6e-11*Params.N_DEG_vec, 'r:');
ylabel('E_f-E_c(k_{||}=0) (meV)');
subplot(212);
plot( Params.N_DEG_vec, k_f_vec./100, 'b', Params.N_DEG_vec, sqrt(Params.N_DEG_vec*2*pi), ':r');
ylabel('k_f (cm^{-1})');
xlabel('N_{2DEG} (cm^{-2})');

%% Optional figure saving

if (Params.save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_structure_profile, [dirPath '\StructureProfile'], 'fig');
        exportfig(h_structure_profile, [dirPath '\StructureProfile'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
    end
end