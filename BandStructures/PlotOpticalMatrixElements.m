function [h_matrix_elements_TE, h_matrix_elements_TM]  = PlotOpticalMatrixElements(k_vec, mu_TE_mat, mu_TM_mat, x, T, L_z, well_material, save)

global Consts;

%% Figures init

h_matrix_elements_TE = figure('Name','Optical Matrix Elements - TE');
h_matrix_elements_TM = figure('Name','Optical Matrix Elements - TM');

%% Figures pre-processing

%E_g = GetMaterialBandGap(well_material,T);

%% Plotting

max_index = ceil(length(k_vec)/2);
[vb_len, cb_len, k_len] = size(mu_TE_mat);

figure(h_matrix_elements_TE);
for (cb_num=1:cb_len)
    subplot(cb_len,1,cb_num); box on; hold on;
    plot(k_vec/100, squeeze(mu_TE_mat(:,cb_num,:)));
    for (vb_num=1:vb_len)
        [m,m_i] = max(squeeze(mu_TE_mat(vb_num,cb_num,:)));
        text(k_vec(m_i), m, ['cb' num2str(cb_num) '-vb' num2str(vb_num)]);
    end
    title(['Optical Matrix elements (TE) - CB' num2str(cb_num)]);
    xlabel('k_{||} (cm^{-1})');
    ylabel('|\mu_{k}^{TE}|^2/|\mu_{k_{0}}|^2');
end

figure(h_matrix_elements_TM);
for (cb_num=1:cb_len)
    subplot(cb_len,1,cb_num); box on; hold on;
    plot(k_vec/100, squeeze(mu_TM_mat(:,cb_num,:)));
    for (vb_num=1:vb_len)
        [m,m_i] = max(squeeze(mu_TM_mat(vb_num,cb_num,:)));
        text(k_vec(m_i), m, ['cb' num2str(cb_num) '-vb' num2str(vb_num)]);
    end
    title(['Optical Matrix elements (TM) - CB' num2str(cb_num)]);
    xlabel('k_{||} (cm^{-1})');
    ylabel('|\mu_{k}^{TM}|^2/|\mu_{k_{0}}|^2');
end

drawnow;

%% Optional figure saving

if (save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_matrix_elements_TE, [dirPath '\OpticalMatrixElement_TE'], 'fig');
        exportfig(h_matrix_elements_TE, [dirPath '\OpticalMatrixElement_TE'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_matrix_elements_TM, [dirPath '\OpticalMatrixElement_TM'], 'fig');
        exportfig(h_matrix_elements_TM, [dirPath '\OpticalMatrixElement_TM'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
    end
end