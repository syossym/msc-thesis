function PlotAbsorptionSpEmission(Structure, Bands, QWParams, MCParams, Params)

global Consts;

[h,w] = size(QWParams);
if (w==1)
    dim = '(:,1)';
    index = '{ii,1}';
elseif(h==1)
    dim = '(1,:)';
    index = '{1,ii}';
end

E_exc = Params.E_exc;
E_g = Structure.QuantumStruct.ActiveLayers{1,1}.E_g;

%% Figures init
scrsz = get(0,'ScreenSize');
pos = [1 0 3*scrsz(3) 1.5*scrsz(4)];

%h_Form_Factors = figure('Name', 'Form Factors');
%h_E_CH = figure('Name', 'Coulomb Hole Self Energy');
%h_E_SX = figure('Name', 'Screened Exchange Shift');
h_Delta_Bandgap = figure('Name', 'Bandgap Corrections');
h_Delta_Bandgap_Total = figure('Name', 'Total Bandgap Corrections');
h_Bandgaps = figure('Name', 'Bandgap');
% h_TE_FCT =  figure('Name', 'Absorption - FCT - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_TE_HF =  figure('Name', 'Absorption - HF - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_TM_FCT =  figure('Name', 'Absorption - FCT - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_TM_HF =  figure('Name', 'Absorption - HF - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_Sp_1 = figure('Name', 'Sp. Emission', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TE_Waterfall = figure('Name', 'Absorption - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TM_Waterfall = figure('Name', 'Absorption - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_Rsp_Waterfall = figure('Name', 'Sp. Emission', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TE_Rsp_Waterfall = figure('Name', 'Sp. Emission vs. Absorption (TE)', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TM_Rsp_Waterfall = figure('Name', 'Sp. Emission vs. Absorption (TM)', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_delta_n_TE_FCT = figure('Name', 'Delta n - FCT - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_delta_n_TM_FCT = figure('Name', 'Delta n - FCT - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_delta_n_TE_HF = figure('Name', 'Delta n - HF - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
% h_delta_n_TM_HF = figure('Name', 'Delta n - HF - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TE_Waterfall = figure('Name', 'Delta n (TE)', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TM_Waterfall = figure('Name', 'Delta n (TM)', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);

%% Variable pre-processing

smooth_param = 3;
norm_alpha = 1e4;
norm_d_n = 1e-4;
norm_rsp = 1e47;

format = '%5.2e';
norm_alpha_str = strrep(num2str(norm_alpha, format), 'e+0', 'x10^{');
norm_alpha_str = strrep(norm_alpha_str, 'e-0', 'x10^{-');
norm_alpha_str = [norm_alpha_str, '}'];
norm_d_n_str = strrep(num2str(norm_d_n, format), 'e+0', 'x10^{');
norm_d_n_str = strrep(norm_d_n_str, 'e-0', 'x10^{-');
norm_d_n_str = [norm_d_n_str, '}'];
norm_rsp_str = strrep(num2str(norm_rsp, format), 'e+0', 'x10^{');
norm_rsp_str = strrep(norm_rsp_str, 'e-0', 'x10^{-');
norm_rsp_str = [norm_rsp_str, '}'];

try
    for (ii=1:eval(['length(QWParams' dim ' )']))
        eval(['current_vec = QWParams' index ';']);
        alpha_TE_FCT{ii} = smooth(current_vec.alpha.TE.FCT./norm_alpha, smooth_param);
        alpha_TM_FCT{ii} = smooth(current_vec.alpha.TM.FCT./norm_alpha, smooth_param);
        alpha_TE_HF{ii} = smooth(current_vec.alpha.TE.HF./norm_alpha, smooth_param);
        alpha_TM_HF{ii} = smooth(current_vec.alpha.TM.HF./norm_alpha, smooth_param);
        alpha_TE_FCT_Waterfall(ii,:) = smooth(current_vec.alpha.TE.FCT, smooth_param);
        alpha_TE_HF_Waterfall(ii,:) = smooth(current_vec.alpha.TE.HF, smooth_param);
        alpha_TM_FCT_Waterfall(ii,:) = smooth(current_vec.alpha.TM.FCT, smooth_param);
        alpha_TM_HF_Waterfall(ii,:) = smooth(current_vec.alpha.TM.HF, smooth_param);
        
        Rsp_FCT{ii} = smooth(current_vec.Rsp.FCT, smooth_param);
        Rsp_HF{ii} = smooth(current_vec.Rsp.HF, smooth_param);
        
        %Rsp_FCT{ii}(Rsp_FCT{ii}<0) = 0;
        %Rsp_HF{ii}(Rsp_HF{ii}<0) = 0;
        %temp_p = Rsp_HF{ii};
        %temp_m = Rsp_HF{ii};
        %temp_p(temp_p<0) = 0; temp_m(temp_m>0) = 0;
        %Rsp_HF{ii} = temp_p-temp_m;
        %Rsp_FCT{ii} = abs(Rsp_FCT{ii});
        Rsp_HF{ii} = abs(Rsp_HF{ii});
        Rsp_FCT{ii} = Rsp_FCT{ii}./norm_rsp;
        Rsp_HF{ii} = Rsp_HF{ii}./norm_rsp;
        Rsp_FCT_Waterfall(ii,:) = Rsp_FCT{ii};
        Rsp_HF_Waterfall(ii,:) = Rsp_HF{ii};
        delta_n_TE_FCT_vec{ii} = smooth(current_vec.delta_n.TE.FCT./norm_d_n, smooth_param);
        delta_n_TM_FCT_vec{ii} = smooth(current_vec.delta_n.TM.FCT./norm_d_n, smooth_param);
        delta_n_TE_FCT_Waterfall(ii,:) = smooth(current_vec.delta_n.TE.FCT, smooth_param);
        delta_n_TM_FCT_Waterfall(ii,:) = smooth(current_vec.delta_n.TM.FCT, smooth_param);
        %delta_n_TE_HF_vec{ii} = abs(delta_n_TE_HF_vec{ii});
        %delta_n_TM_HF_vec{ii} = abs(delta_n_TM_HF_vec{ii});
        delta_n_TE_HF_vec{ii} = smooth(current_vec.delta_n.TE.HF./norm_d_n, smooth_param);
        delta_n_TM_HF_vec{ii} = smooth(current_vec.delta_n.TM.HF./norm_d_n, smooth_param);
        delta_n_TE_HF_Waterfall(ii,:) = smooth(current_vec.delta_n.TE.HF, smooth_param);
        delta_n_TM_HF_Waterfall(ii,:) = smooth(current_vec.delta_n.TM.HF, smooth_param);
    end
catch exc0
    disp(exc0.message);
    %conitnue
end

%% Plotting defintions

num_plots = length(alpha_TE_FCT_Waterfall(:,1));
n = 4;
num_subplot_rows = ceil(num_plots/n);
%E_grid = (E_exc/Consts.e_0-E_g)*1e3;     % (E-E_g) [meV]
E_grid = E_exc/Consts.e_0;                % E (eV)
E_11 = Bands{1,1}.Cond.E_0(1) - Bands{1,1}.Valence.E_0(1);
E_s = 0;     % E_g                        % shift energy for the grid
plot_grid_min = 1.52; % min(E_grid);     % (eV)
plot_grid_max = 1.55; % max(E_grid);       % (eV)

%% Plotting

try
    for (ii=1:num_plots)
        delta_E_CH_k_0(1,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_CH{1}(1)/Consts.e_0;
        delta_E_CH_k_0(2,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_CH{2}(1)/Consts.e_0;
        delta_E_CH_k_0(3,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_CH{3}(1)/Consts.e_0;
        delta_E_SX_k_0(1,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_SX{1}(1)/Consts.e_0;
        delta_E_SX_k_0(2,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_SX{2}(1)/Consts.e_0;
        delta_E_SX_k_0(3,ii) = 1e3*QWParams{ii}.CalculatedParams.D_E_SX{3}(1)/Consts.e_0;
    end
    
    figure(h_Delta_Bandgap);
    subplot(211); box on;
    semilogx(Params.N_DEG_vec, delta_E_CH_k_0.');
    ylabel('\DeltaE_{CH,n}(k_{||}=0) (meV)'); title('(a)'); set(gca, 'XTickLabel', '');
    legend('n=1','n=2','n=3');
    subplot(212); box on;
    semilogx(Params.N_DEG_vec, delta_E_SX_k_0.');
    ylabel('\DeltaE_{CH,nm}(k_{||}=0) (meV)'); title('(b)');
    xlabel('N_{2DEG} (cm^{-2})');
    legend('n,m=1,1','n,m=1,2','n,m=1,3');
    
    figure(h_Delta_Bandgap_Total);
    semilogx(Params.N_DEG_vec,(delta_E_CH_k_0+delta_E_SX_k_0).');
    ylabel('Bandgap Correction (k_{||}=0) (meV)'); xlabel('N_{2DEG} (cm^{-2})');
    legend('n=1','n=2','n=3');
    
    k_vec = interp(Bands{1}.k_t_vec(Params.k_indices), Params.k_scale);
    figure(h_Bandgaps);
    for (ii=1:num_plots)
        subplot(311); hold on; box on;
        plot(k_vec/100, QWParams{ii}.CalculatedParams.w_mn{1}*Consts.hbar/Consts.e_0, 'b', k_vec/100, QWParams{ii}.CalculatedParams.w_mn_tag{1}*Consts.hbar/Consts.e_0, 'r', 'LineWidth', 0.1);
        ylabel('\hbar\omega_{n,m}, \hbar\tilde{\omega_{n,m}} (eV)'); title('(a)'); set(gca, 'XTickLabel', '');
        subplot(312); hold on; box on;
        plot(k_vec/100, QWParams{ii}.CalculatedParams.w_mn{2}*Consts.hbar/Consts.e_0, 'b', k_vec/100, QWParams{ii}.CalculatedParams.w_mn_tag{2}*Consts.hbar/Consts.e_0, 'r', 'LineWidth', 0.1);
        title('(b)'); set(gca, 'XTickLabel', '');
        subplot(313); hold on; box on;
        plot(k_vec/100, QWParams{ii}.CalculatedParams.w_mn{3}*Consts.hbar/Consts.e_0, 'b', k_vec/100, QWParams{ii}.CalculatedParams.w_mn_tag{3}*Consts.hbar/Consts.e_0, 'r', 'LineWidth', 0.1);
        xlabel('k_{||} (cm^{-1})'); title('(c)');
    end
catch exc
    
end

% try
%     for (ii=1:num_plots)
%         figure(h_TE_FCT);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (alpha_TE_FCT{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TE_FCT{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(alpha_TE_FCT{ii}), max(alpha_TE_FCT{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\alpha_{TE}^{FCT} (' norm_alpha_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         [maxs, mins] = peakdet(alpha_TE_HF{ii}, 0.1, E_exc/Consts.e_0);
%         
%         figure(h_TE_HF);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (alpha_TE_HF{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TE_HF{ii}, 'k-');   % E_g0
%         for (mx=1:length(maxs(:,1)))
%             text(maxs(mx,1), maxs(mx,2), [num2str(maxs(mx,1)) 'eV']);
%         end
%         axis([plot_grid_min, plot_grid_max, min(alpha_TE_HF{ii}), max(alpha_TE_HF{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\alpha_{TE}^{HF} (' norm_alpha_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         figure(h_TM_FCT);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (alpha_TM_FCT{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TM_FCT{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(alpha_TM_FCT{ii}), max(alpha_TM_FCT{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\alpha_{TM}^{FCT} (' norm_alpha_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         [maxs, mins] = peakdet(alpha_TM_HF{ii}, 0.1, E_exc/Consts.e_0);
%         
%         figure(h_TM_HF);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (alpha_TM_HF{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TM_HF{ii}, 'k-');   % E_g0
%         for (mx=1:length(maxs(:,1)))
%             text(maxs(mx,1), maxs(mx,2), [num2str(maxs(mx,1)) 'eV']);
%         end
%         axis([plot_grid_min, plot_grid_max, min(alpha_TM_HF{ii}), max(alpha_TM_HF{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\alpha_{TM}^{HF} (' norm_alpha_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         figure(h_delta_n_TE_FCT);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (delta_n_TE_FCT_vec{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TE_FCT_vec{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(delta_n_TE_FCT_vec{ii}), max(delta_n_TE_FCT_vec{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\deltan_{TE}^{FCT} (' norm_d_n_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     %xlabel('E-E_g [meV]');
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         figure(h_delta_n_TM_FCT);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (delta_n_TM_FCT_vec{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TM_FCT_vec{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(delta_n_TM_FCT_vec{ii}), max(delta_n_TM_FCT_vec{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\deltan_{TM}^{FCT} (' norm_d_n_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         figure(h_delta_n_TE_HF);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (delta_n_TE_HF_vec{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TE_HF_vec{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(delta_n_TE_HF_vec{ii}), max(delta_n_TE_HF_vec{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\deltan_{TE}^{HF} (' norm_d_n_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         figure(h_delta_n_TM_HF);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, (delta_n_TM_HF_vec{ii}));
%         plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TM_HF_vec{ii}, 'k-');   % E_g0
%         axis([plot_grid_min, plot_grid_max, min(delta_n_TM_HF_vec{ii}), max(delta_n_TM_HF_vec{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['\deltan_{TM}^{HF} (' norm_d_n_str ' cm^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end
% 
% try
%     for (ii=1:num_plots)
%         [maxs_HF, mins_HF] = peakdet(Rsp_HF{ii}, 0.1, E_exc/Consts.e_0);
%         
%         figure(h_Sp_1);
%         subplot(num_plots,1,ii); hold on; box on; axis tight;
%         plot(E_grid, Rsp_FCT{ii});
%         plot(E_grid, Rsp_HF{ii}, 'r');
%         plot((E_g-E_s)*ones(1,length(E_grid)), Rsp_HF{ii}, 'k-');   % E_g0
%         [max_FCT, max_FCT_index] = max(Rsp_FCT{ii});
%         text(E_grid(max_FCT_index), max_FCT, [num2str(E_grid(max_FCT_index)) 'eV']);
%         for (mx=1:length(maxs_HF(:,1)))
%             text(maxs_HF(mx,1), maxs_HF(mx,2), [num2str(maxs_HF(mx,1)) 'eV']);
%         end
%         axis([plot_grid_min, plot_grid_max, min(Rsp_HF{ii}), max(Rsp_HF{ii})]);
%         if (w==1)
%             y_text = strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         elseif (h==1)
%             y_text = strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{');
%         end
%         y_text = [y_text '}'];
%         ylabel(y_text);
%         h_curr = gca;
%         if (ii==1)
%             title(['R_{sp}^{FCT} (blue), R_{sp}^{HF} (red) (' norm_rsp_str ' s^{-1}m^{-3}J^{-1})']);
%         end
%         if (ii~=num_plots)
%             set(h_curr, 'XTickLabel', '');
%         end
%     end
%     xlabel('E (eV)');
% catch exc4
%     disp(exc4.message);
%     % conitnue
% end

if (w==1)
    y_var.name = 'N_{2DEG}';
    y_var.value = Params.N_DEG_vec;
    y_var.units = 'cm^{-2}';
    y_var.num_type = 'exp';
elseif (h==1)
    y_var.name = '\gamma';
    y_var.value = Params.gamma_vec;
    y_var.units = 'sec^{-1}';
    y_var.num_type = 'exp';
end

grid_range = [plot_grid_min, plot_grid_max];

PlotWaterfall(E_grid, alpha_TE_FCT_Waterfall, h_TE_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, alpha_TE_HF_Waterfall, h_TE_Waterfall, 'red', grid_range);
ylabel('\alpha^{FCT}_{TE} (blue), \alpha^{HF}_{TE} (red) (a.u.)');
xlabel('E (eV)');
%title('Absorption spectra - TE polarization');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, alpha_TM_FCT_Waterfall, h_TM_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, alpha_TM_HF_Waterfall, h_TM_Waterfall, 'red', grid_range);
ylabel('\alpha^{FCT}_{TM} (blue), \alpha^{HF}_{TM} (red) (a.u.)');
xlabel('E (eV)');
%title('Absorption spectra - TM polarization');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, Rsp_FCT_Waterfall, h_Rsp_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, Rsp_HF_Waterfall, h_Rsp_Waterfall, 'red', grid_range);
ylabel('r^{FCT}_{sp} (blue), r^{HF}_{sp} (red) (a.u.)');
xlabel('E (eV)');
%title('Spontanious Emission');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, alpha_TE_HF_Waterfall, h_TE_Rsp_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, Rsp_HF_Waterfall, h_TE_Rsp_Waterfall, 'red', grid_range);
ylabel('\alpha^{HF}_{TE} (blue), r^{HF}_{sp} (red) (a.u.)');
xlabel('E (eV)');
%title('Spontanious Emission vs. Absorption (TE)');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, alpha_TM_HF_Waterfall, h_TM_Rsp_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, Rsp_HF_Waterfall, h_TM_Rsp_Waterfall, 'red', grid_range);
ylabel('\alpha^{HF}_{TM} (blue), r^{HF}_{sp} (red) (a.u.)');
xlabel('E (eV)');
%title('Spontanious Emission vs. Absorption (TM)');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, delta_n_TE_FCT_Waterfall, h_delta_n_TE_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, delta_n_TE_HF_Waterfall, h_delta_n_TE_Waterfall, 'red', grid_range);
ylabel('\deltan^{FCT}_{TE} (blue), \deltan^{HF}_{TE} (red) (a.u.)');
xlabel('E (eV)');
%title('\deltan (TE)');
set(gca, 'YTickLabel', '');

PlotWaterfall(E_grid, delta_n_TM_FCT_Waterfall, h_delta_n_TM_Waterfall, 'blue', grid_range, y_var);
PlotWaterfall(E_grid, delta_n_TM_HF_Waterfall, h_delta_n_TM_Waterfall, 'red', grid_range);
ylabel('\deltan^{FCT}_{TM} (blue), \deltan^{HF}_{TM} (red) (a.u.)');
xlabel('E (eV)');
%title('\deltan (TM)');
set(gca, 'YTickLabel', '');

%% Optional figure saving

if (Params.save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_TE_FCT, [dirPath '\Abs_FCT_TE'], 'fig');
        exportfig(h_TE_FCT, [dirPath '\Abs_FCT_TE'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_TE_HF, [dirPath '\Abs_HF_TE'], 'fig');
        exportfig(h_TE_HF, [dirPath '\Abs_HF_TE'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_TM_FCT, [dirPath '\Abs_FCT_TM'], 'fig');
        exportfig(h_TM_FCT, [dirPath '\Abs_FCT_TM'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_TM_HF, [dirPath '\Abs_HF_TM'], 'fig');
        exportfig(h_TM_HF, [dirPath '\Abs_HF_TM'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_Sp_1, [dirPath '\Sp_All'], 'fig');
        exportfig(h_Sp_1, [dirPath '\Sp_All'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_delta_n_TE_FCT, [dirPath '\Delta_n_FCT_TE'], 'fig');
        exportfig(h_delta_n_TE_FCT, [dirPath '\Delta_n_FCT_TE'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_delta_n_TM_FCT, [dirPath '\Delta_n_FCT_TM'], 'fig');
        exportfig(h_delta_n_TM_FCT, [dirPath '\Delta_n_FCT_TM'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_delta_n_TE_HF, [dirPath '\Delta_n_HF_TE'], 'fig');
        exportfig(h_delta_n_TE_HF, [dirPath '\Delta_n_HF_TE'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
        saveas(h_delta_n_TM_HF, [dirPath '\Delta_n_HF_TM'], 'fig');
        exportfig(h_delta_n_TM_HF, [dirPath '\Delta_n_HF_TM'], 'FontMode', 'scaled', 'FontSize', 0.5, 'color', 'cmyk', 'Format', 'eps');
    end
end