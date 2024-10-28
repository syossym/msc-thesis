function [h_Abs_Sp_FCT, h_Abs_Sp_HF] = PlotAbsorptionSpEmission(E_exc, alpha_TE_FCT, alpha_TM_FCT, alpha_TE_HF, alpha_TM_HF, Rsp_FCT, Rsp_HF, delta_n_TE_FCT_vec, delta_n_TM_FCT_vec, delta_n_TE_HF_vec, delta_n_TM_HF_vec, N_e_vec, T, E_f_vec, E_c, E_v, delta_renorm, E_g, save)

global Consts;

%% Figures init
scrsz = get(0,'ScreenSize');
pos = [1 0 3*scrsz(3) 1.5*scrsz(4)];

h_TE_FCT =  figure('Name', 'Absorption - FCT - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TE_HF =  figure('Name', 'Absorption - HF - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TM_FCT =  figure('Name', 'Absorption - FCT - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TM_HF =  figure('Name', 'Absorption - HF - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_Sp_1 = figure('Name', 'Sp. Emission', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TE_Waterfall = figure('Name', 'Absorption - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_TM_Waterfall = figure('Name', 'Absorption - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_Rsp_Waterfall = figure('Name', 'Sp. Emission', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TE_FCT = figure('Name', 'Delta n - FCT - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TM_FCT = figure('Name', 'Delta n - FCT - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TE_HF = figure('Name', 'Delta n - HF - TE', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);
h_delta_n_TM_HF = figure('Name', 'Delta n - HF - TM', 'PaperPosition',[4 4 13 20],'Position',[0 0 480 500]);


%% Variable pre-processing

norm_alpha = 1e4;
norm_d_n = 1e-4;
norm_rsp = 1e47;

format = '%1.2e';
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
    for (ii=1:length(N_e_vec))
        alpha_TE_FCT{ii} = alpha_TE_FCT{ii}./norm_alpha;
        alpha_TM_FCT{ii} = alpha_TM_FCT{ii}./norm_alpha;
        alpha_TE_HF{ii} = alpha_TE_HF{ii}./norm_alpha;
        alpha_TM_HF{ii} = alpha_TM_HF{ii}./norm_alpha;
        alpha_TE_FCT_Waterfall(ii,:) = alpha_TE_FCT{ii};
        alpha_TE_HF_Waterfall(ii,:) = alpha_TE_HF{ii};
        alpha_TM_FCT_Waterfall(ii,:) = alpha_TM_FCT{ii};
        alpha_TM_HF_Waterfall(ii,:) = alpha_TM_HF{ii};
        Rsp_FCT{ii}(Rsp_FCT{ii}<0) = 0;
        Rsp_HF{ii}(Rsp_HF{ii}<0) = 0;
        %temp_p = Rsp_HF{ii};
        %temp_m = Rsp_HF{ii};
        %temp_p(temp_p<0) = 0; temp_m(temp_m>0) = 0;
        %Rsp_HF{ii} = temp_p-temp_m;
        %Rsp_FCT{ii} = abs(Rsp_FCT{ii});
        %Rsp_HF{ii} = abs(Rsp_HF{ii});
        Rsp_FCT{ii} = Rsp_FCT{ii}./norm_rsp;
        Rsp_HF{ii} = Rsp_HF{ii}./norm_rsp;
        Rsp_FCT_Waterfall(ii,:) = Rsp_FCT{ii};
        Rsp_HF_Waterfall(ii,:) = Rsp_HF{ii};
        delta_n_TE_FCT_vec{ii} = delta_n_TE_FCT_vec{ii}./norm_d_n;
        delta_n_TM_FCT_vec{ii} = delta_n_TM_FCT_vec{ii}./norm_d_n;
        %delta_n_TE_HF_vec{ii} = abs(delta_n_TE_HF_vec{ii});
        %delta_n_TM_HF_vec{ii} = abs(delta_n_TM_HF_vec{ii});
        delta_n_TE_HF_vec{ii} = delta_n_TE_HF_vec{ii}./norm_d_n;
        delta_n_TM_HF_vec{ii} = delta_n_TM_HF_vec{ii}./norm_d_n;
    end
catch exc0
    disp(exc0.message);
    %conitnue
end


%% Plotting defintions

num_cons = length(alpha_TE_FCT_Waterfall(:,1));
n = 4;
num_subplot_rows = ceil(num_cons/n);
%E_grid = (E_exc/Consts.e_0-E_g)*1e3;     % (E-E_g) [meV]
E_grid = E_exc/Consts.e_0;                % E [eV]
E_11 = E_c{1}(1) + E_v{1}(1);
E_s = 0;     % E_g                        % shift energy for the grid
plot_grid_min = 1.515; % min(E_grid);     % [eV]
plot_grid_max = 1.54; %max(E_grid);       % [eV]


%% Plotting

try
    for (ii=1:num_cons)
        figure(h_TE_FCT);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (alpha_TE_FCT{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TE_FCT{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_FCT{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_FCT{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_FCT{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_FCT{ii}, 'g:');        % min(E_c_1) renorm
        axis([plot_grid_min, plot_grid_max, min(alpha_TE_FCT{ii}), max(alpha_TE_FCT{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\alpha_{TE}^{FCT} (' norm_alpha_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        [maxs, mins] = peakdet(alpha_TE_HF{ii}, 0.1, E_exc/Consts.e_0);
        
        figure(h_TE_HF);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (alpha_TE_HF{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TE_HF{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_HF{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_HF{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_HF{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TE_HF{ii}, 'g:');        % min(E_c_1)
        for (mx=1:length(maxs(:,1)))
            text(maxs(mx,1), maxs(mx,2), [num2str(maxs(mx,1)) 'eV']);
        end
        axis([plot_grid_min, plot_grid_max, min(alpha_TE_HF{ii}), max(alpha_TE_HF{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\alpha_{TE}^{HF} (' norm_alpha_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        figure(h_TM_FCT);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        plot(E_grid, (alpha_TM_FCT{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TM_FCT{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_FCT{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_FCT{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_FCT{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_FCT{ii}, 'g:');        % min(E_c_1)
        axis([plot_grid_min, plot_grid_max, min(alpha_TM_FCT{ii}), max(alpha_TM_FCT{ii})]);
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\alpha_{TM}^{FCT} (' norm_alpha_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        [maxs, mins] = peakdet(alpha_TM_HF{ii}, 0.1, E_exc/Consts.e_0);
        
        figure(h_TM_HF);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        plot(E_grid, (alpha_TM_HF{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), alpha_TM_HF{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_HF{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_HF{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_HF{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), alpha_TM_HF{ii}, 'g:');        % min(E_c_1)
        for (mx=1:length(maxs(:,1)))
            text(maxs(mx,1), maxs(mx,2), [num2str(maxs(mx,1)) 'eV']);
        end
        axis([plot_grid_min, plot_grid_max, min(alpha_TM_HF{ii}), max(alpha_TM_HF{ii})]);
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\alpha_{TM}^{HF} (' norm_alpha_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        figure(h_delta_n_TE_FCT);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (delta_n_TE_FCT_vec{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TE_FCT_vec{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_FCT_vec{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_FCT_vec{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_FCT_vec{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_FCT_vec{ii}, 'g:');        % min(E_c_1) renorm
        axis([plot_grid_min, plot_grid_max, min(delta_n_TE_FCT_vec{ii}), max(delta_n_TE_FCT_vec{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\deltan_{TE}^{FCT} (' norm_d_n_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        figure(h_delta_n_TM_FCT);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (delta_n_TM_FCT_vec{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TM_FCT_vec{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_FCT_vec{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_FCT_vec{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_FCT_vec{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_FCT_vec{ii}, 'g:');        % min(E_c_1) renorm
        axis([plot_grid_min, plot_grid_max, min(delta_n_TM_FCT_vec{ii}), max(delta_n_TM_FCT_vec{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\deltan_{TM}^{FCT} (' norm_d_n_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        figure(h_delta_n_TE_HF);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (delta_n_TE_HF_vec{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TE_HF_vec{ii}, 'k-');   % E_g0
        %plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_HF_vec{ii}, 'r-');      % E_f
        %plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_HF_vec{ii}, 'r:');      % E_f renorm
        %plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_HF_vec{ii}, 'g-');        % min(E_c_1)
        %plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TE_HF_vec{ii}, 'g:');        % min(E_c_1) renorm
        axis([plot_grid_min, plot_grid_max, min(delta_n_TE_HF_vec{ii}), max(delta_n_TE_HF_vec{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\deltan_{TE}^{HF} (' norm_d_n_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        figure(h_delta_n_TM_HF);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        %[AX,H1,H2] = plotyy(E_grid, alpha_TE_FCT{ii}, E_grid, alpha_TE_HF{ii});
        plot(E_grid, (delta_n_TM_HF_vec{ii}));
        plot((E_g-E_s)*ones(1,length(E_grid)), delta_n_TM_HF_vec{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_HF_vec{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_HF_vec{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_HF_vec{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), delta_n_TM_HF_vec{ii}, 'g:');        % min(E_c_1) renorm
        axis([plot_grid_min, plot_grid_max, min(delta_n_TM_HF_vec{ii}), max(delta_n_TM_HF_vec{ii})]);
        %set(get(AX(1), 'Ylabel'), 'String', [num2str(N_e_vec(ii), format) 'cm^{-2}']);
        %set(H1,'LineStyle','-', 'Color', 'r')
        %set(H2,'LineStyle','-', 'Color', 'b')
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['\deltan_{TM}^{HF} (' norm_d_n_str ' cm^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

try
    for (ii=1:num_cons)
        [maxs_HF, mins_HF] = peakdet(Rsp_HF{ii}, 0.1, E_exc/Consts.e_0);
        
        figure(h_Sp_1);
        subplot(num_cons,1,ii); hold on; box on; axis tight;
        plot(E_grid, Rsp_FCT{ii});
        plot(E_grid, Rsp_HF{ii}, 'r');
        plot((E_g-E_s)*ones(1,length(E_grid)), Rsp_HF{ii}, 'k-');   % E_g0
        %         plot((E_f_vec(ii)*ones(1,length(E_grid))/Consts.e_0-E_s), Rsp_HF{ii}, 'r-');      % E_f
        %         plot(((E_f_vec(ii)+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), Rsp_HF{ii}, 'r:');      % E_f renorm
        %         plot((E_11*ones(1,length(E_grid))/Consts.e_0-E_s), Rsp_HF{ii}, 'g-');        % min(E_c_1)
        %         plot(((E_11+delta_renorm(ii))*ones(1,length(E_grid))/Consts.e_0-E_s), Rsp_HF{ii}, 'g:');        % min(E_c_1)
        [max_FCT, max_FCT_index] = max(Rsp_FCT{ii});
        text(E_grid(max_FCT_index), max_FCT, [num2str(E_grid(max_FCT_index)) 'eV']);
        for (mx=1:length(maxs_HF(:,1)))
            text(maxs_HF(mx,1), maxs_HF(mx,2), [num2str(maxs_HF(mx,1)) 'eV']);
        end
        axis([plot_grid_min, plot_grid_max, min(Rsp_HF{ii}), max(Rsp_HF{ii})]);
        con_text = strrep([num2str(N_e_vec(ii), format)], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        if (ii==1)
            title(['R_{sp}^{FCT} (blue), R_{sp}^{HF} (red) (' norm_rsp_str ' s^{-1}m^{-3}J^{-1})']);
        end
        if (ii~=num_cons)
            set(h_curr, 'XTickLabel', '');
        end
    end
    %xlabel('E-E_g [meV]');
    xlabel('E [eV]');
catch exc4
    disp(exc4.message);
    % conitnue
end

figure(h_TE_Waterfall);
subplot(211);
plot(E_grid, alpha_TE_FCT_Waterfall);
ylabel('\alpha^{FCT}_{TE} [cm^{-1}]');
title('Absorption spectra - TE polarization');
axis([plot_grid_min, plot_grid_max, min(min(alpha_TE_FCT_Waterfall)), max(max(alpha_TE_FCT_Waterfall))]);
subplot(212); hold on; box on;
plot(E_grid, alpha_TE_HF_Waterfall);
%xlabel('E-E_g [meV]');
xlabel('E [eV]');
ylabel('\alpha^{HF}_{TE} [cm^{-1}]');
axis([plot_grid_min, plot_grid_max, min(min(alpha_TE_HF_Waterfall)), max(max(alpha_TE_HF_Waterfall))]);

figure(h_TM_Waterfall);
subplot(211);
plot(E_grid, alpha_TM_FCT_Waterfall);
ylabel('\alpha^{FCT}_{TM} [cm^{-1}]');
title('Absorption spectra - TM polarization');
axis([plot_grid_min, plot_grid_max, min(min(alpha_TM_FCT_Waterfall)), max(max(alpha_TM_FCT_Waterfall))]);
subplot(212); hold on; box on;
plot(E_grid, alpha_TM_HF_Waterfall);
%xlabel('E-E_g [meV]');
xlabel('E [eV]');
ylabel('\alpha^{HF}_{TM} [cm^{-1}]');
axis([plot_grid_min, plot_grid_max, min(min(alpha_TM_HF_Waterfall)), max(max(alpha_TM_HF_Waterfall))]);

figure(h_Rsp_Waterfall);
subplot(211);
plot(E_grid, Rsp_FCT_Waterfall);
ylabel('R^{FCT}_{sp} [a.u.]');
title('Spontanious Emission');
axis([plot_grid_min, plot_grid_max, min(min(Rsp_FCT_Waterfall)), max(max(Rsp_FCT_Waterfall))]);
subplot(212); hold on; box on;
plot(E_grid, Rsp_HF_Waterfall);
%xlabel('E-E_g [meV]');
xlabel('E [eV]');
ylabel('R^{HF}_{sp} [a.u.]');
axis([plot_grid_min, plot_grid_max, min(min(Rsp_HF_Waterfall)), max(max(Rsp_HF_Waterfall))]);


%% Optional figure saving

if (save)
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