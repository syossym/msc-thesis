function PlotOpticalSimulationParams(Structure, Bands, QWParams, MCParams, Params)
%
% This function plots the optical simulation results.
%
%   Input: 'Structure' - the simulated structure.
%          'Bands' - the quntum band structure.
%          'QWParams' - the results of the optical parameters calculation.
%          'MCParams' - the results of the MC structure reflection calculation.
%          'Params' - simulaion parameters.
%
%   Output: -
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

y_var.num_type = 'exp';

h_RefractiveIndexProfile = figure('Name', 'Refractive Index Profile');

figure(h_RefractiveIndexProfile);
plot(MCParams{1,1}.z_grid/1e-10, MCParams{1,1}.n_profile); hold on;
axis([min(MCParams{1,1}.z_grid)/1e-10, max(MCParams{1,1}.z_grid)/1e-10, 0 , max(MCParams{1,1}.n_profile)+0.2]);
xlabel('z [A]'); ylabel('n');

if (length(Params.N_DEG_vec) > 1 && length(Params.gamma_vec) == 1)
    sim_type = 1;
elseif (length(Params.N_DEG_vec) == 1 && length(Params.gamma_vec) > 1)
    sim_type = 2;
elseif (length(Params.N_DEG_vec) == 1 && length(Params.gamma_vec) == 1)
    sim_type = 3;
end

colors_vec = ['b', 'r', 'g', 'k'];
reflection_resolution = 1;

h_Reflection_TE = figure('Name', 'Reflection (TE)');
h_Reflection_TM = figure('Name', 'Reflection (TM)');

h_AntiCrossing_TE = figure('Name', 'Anticrossing Diagrams (TE)');
h_AntiCrossing_TM = figure('Name', 'Anticrossing Diagrams (TM)');
h_AntiCrossing_TE_TM = figure('Name', 'Anticrossing Diagrams (TE & TM)');
h_AntiCrossing_TE_Waterfall = figure('Name', 'Anticrossing Diagrams (TE)');
h_AntiCrossing_TM_Waterfall = figure('Name', 'Anticrossing Diagrams (TM)');
h_AntiCrossing_Waterfall = figure('Name', 'Anticrossing Diagrams');
h_Susceptability_TE = figure('Name', 'QW Susceptibility (TE)');
h_Susceptability_TM = figure('Name', 'QW Susceptibility (TM)');
h_n_TE = figure('Name', 'QW Refractive Index (TE)');
h_n_TM = figure('Name', 'QW Refractive Index (TM)');
h_Susceptability_TE_Waterfall = figure('Name', 'QW Refractive Index (TE)');
h_Susceptability_TM_Waterfall = figure('Name', 'QW Refractive Index (TM)');

switch (sim_type)
    
    case 1,    % multiple DEG concentrations and single gamma value
        
        for (ii=1:length(MCParams(:,1)))
            
            y_var.name = '\delta';
            y_var.value = nan; %  MCParams{1,1}.Detuning.delta_vec;
            y_var.units = '';
            grid_range = [1.5, 1.55];
            
            arranged_curves_TE = ArrangeAntiCrossingData(MCParams{ii,1}.Detuning.E_min_MC,MCParams{ii,1}.Detuning.E_min_MCQW_TE);
            arranged_curves_TM = ArrangeAntiCrossingData(MCParams{ii,1}.Detuning.E_min_MC,MCParams{ii,1}.Detuning.E_min_MCQW_TM);
%             for (jj=1:length(arranged_curves_TE))
%                 if (~isempty(arranged_curves_TE{jj}))
%                     figure(h_AntiCrossing_TE);
%                     subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
%                     plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MCQW, 'b.', 'MarkerSize', 2);
%                     plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MC, 'r');
%                     %xlabel('E_{MC} (eV)'); ylabel('E (eV)');
%                     title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
%                     axis([grid_range grid_range]);
%                     figure(h_AntiCrossing_TE_TM);
%                     subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
%                     plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MCQW, 'b.', 'MarkerSize', 2);
%                 end
%             end
%             for (jj=1:length(arranged_curves_TM))
%                 if (~isempty(arranged_curves_TM{jj}))
%                     figure(h_AntiCrossing_TM);
%                     subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
%                     plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MCQW, 'b.', 'MarkerSize', 2);
%                     plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MC, 'r');
%                     %xlabel('E_{MC} (eV)'); ylabel('E (eV)');
%                     title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);                  
%                     axis([grid_range grid_range]);
%                     figure(h_AntiCrossing_TE_TM);
%                     subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
%                     plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MCQW, 'g.', 'MarkerSize', 2);
%                 end
%             end
            
            figure(h_AntiCrossing_TE_TM);
            subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            plot(arranged_curves_TM{1}.E_MC, arranged_curves_TM{1}.E_MC, 'r');
            title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
            axis([grid_range grid_range]);
                       
            reflection_vector = 1:round(length(MCParams{ii,1}.Detuning.r_MCQW_TE(:,1))/reflection_resolution):length(MCParams{ii,1}.Detuning.r_MCQW_TE(:,1));
            
            figure(h_AntiCrossing_TE_Waterfall);
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            PlotWaterfall(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, abs(MCParams{ii,1}.Detuning.r_MCQW_TE(reflection_vector,:)), h, 'b', grid_range, y_var, 5);
            title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_AntiCrossing_TM_Waterfall);
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            PlotWaterfall(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, MCParams{ii,1}.Detuning.r_MCQW_TM(reflection_vector,:), h, 'b', grid_range, y_var, 5);
            title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_AntiCrossing_Waterfall); hold on;
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            PlotWaterfall(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, MCParams{ii,1}.Detuning.r_MCQW_TE(reflection_vector,:), h, 'b', grid_range, y_var, 5);
            PlotWaterfall(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, MCParams{ii,1}.Detuning.r_MCQW_TM(reflection_vector,:), h, 'r', grid_range, y_var, 5);
            title(['N_{2DEG}=' strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_n_TE);
            subplot(length(MCParams(:,1)),1,ii);
            [AX,H1,H2] = plotyy(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, smooth(real(MCParams{ii,1}.Detuning.n_QW_TE), 30), MCParams{ii,1}.Detuning.E_vec/Consts.e_0, smooth(imag(MCParams{ii,1}.Detuning.n_QW_TE), 30));
            %set(get(AX(2), 'Ylabel'), 'String', '\Im[n_{QW}]');
            set(H1,'LineStyle','-', 'Color', 'b');
            set(H2,'LineStyle','-', 'Color', 'r');
            set(AX(1), 'XLim', [min(MCParams{ii,1}.Detuning.E_vec/Consts.e_0), max(MCParams{ii,1}.Detuning.E_vec/Consts.e_0)]); %, 'YTick', [3.5, 5]);
            set(AX(2), 'XLim', [min(MCParams{ii,1}.Detuning.E_vec/Consts.e_0), max(MCParams{ii,1}.Detuning.E_vec/Consts.e_0)]); %, 'YTick', [0, 1]);
            %set(AX(1), 'YLim', [3.5, max(real(MCParams{1,1}.Detuning.n_QW_TE))]);
            %set(AX(2), 'YLim', [0, max(imag(MCParams{1,1}.Detuning.n_QW_TE))]);
            set(AX(1), 'YColor', 'b');
            set(AX(2), 'YColor', 'r');
            con_text = strrep([num2str(MCParams{ii,1}.N_DEG, '%1.0e')], 'e+0', 'x10^{');
            con_text = [con_text '} cm^{-2}'];
            ylabel(con_text);
            h_curr = gca;
            if (ii==1)
                title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
                title(title_text);
            end
            if (ii~=length(MCParams(:,1)))
                set(AX(1), 'XTickLabel', '');
                set(AX(2), 'XTickLabel', '');
            end
            
            figure(h_n_TM);
            subplot(length(MCParams(:,1)),1,ii);
            [AX,H1,H2] = plotyy(MCParams{ii,1}.Detuning.E_vec/Consts.e_0, smooth(real(MCParams{ii,1}.Detuning.n_QW_TM), 30), MCParams{ii,1}.Detuning.E_vec/Consts.e_0, smooth(imag(MCParams{ii,1}.Detuning.n_QW_TM), 30));
            %set(get(AX(2), 'Ylabel'), 'String', '\Im[n_{QW}]');
            axis auto; axis tight;
            set(H1,'LineStyle','-', 'Color', 'b');
            set(H2,'LineStyle','-', 'Color', 'r');
            set(AX(1), 'XLim', [min(MCParams{ii,1}.Detuning.E_vec/Consts.e_0), max(MCParams{ii,1}.Detuning.E_vec/Consts.e_0)]); %, 'YTick', [3.5, 5]);
            set(AX(2), 'XLim', [min(MCParams{ii,1}.Detuning.E_vec/Consts.e_0), max(MCParams{ii,1}.Detuning.E_vec/Consts.e_0)]); %, 'YTick', [0, 1]);
            %set(AX(1), 'YLim', [3.5, max(real(MCParams{1,1}.Detuning.n_QW_TM))]);
            %set(AX(2), 'YLim', [0, max(imag(MCParams{1,1}.Detuning.n_QW_TM))]);
            set(AX(1), 'YColor', 'b');
            set(AX(2), 'YColor', 'r');
            con_text = strrep([num2str(MCParams{ii,1}.N_DEG, '%1.0e')], 'e+0', 'x10^{');
            con_text = [con_text '} cm^{-2}'];
            ylabel(con_text);
            h_curr = gca;
            if (ii==1)
                title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
                title(title_text);
            end
            if (ii~=length(MCParams(:,1)))
                set(AX(1), 'XTickLabel', '');
                set(AX(2), 'XTickLabel', '');
            end
            
            n_QW_TE_Waterfall(ii,:) = MCParams{ii,1}.Detuning.n_QW_TE;
            n_QW_TM_Waterfall(ii,:) = MCParams{ii,1}.Detuning.n_QW_TM;
                       
        end
        
        figure(h_Susceptability_TE);
        xlabel('E (eV)');
        figure(h_Susceptability_TM);
        xlabel('E (eV)');
        
        y_var.name = 'N_{2GED}';
        y_var.value = Params.N_DEG_vec;
        y_var.units = 'cm^{-2}';
        
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, real(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'b', grid_range, y_var);
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['Refractive Index (TE): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');
        
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, real(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'b', grid_range, y_var);
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['Refractive Index (TM): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');

    case 2,    % single DEG concentration and multiple gamma values
        
        for (ii=1:length(MCParams(1,:)))
            
            y_var.name = '\delta';
            y_var.value = nan;%MCParams{1,1}.Detuning.delta_vec;
            y_var.units = '';
            grid_range = [1.52, 1.54];
            
            arranged_curves_TE = ArrangeAntiCrossingData(MCParams{1,ii}.Detuning.E_min_MC,MCParams{1,ii}.Detuning.E_min_MCQW_TE);
            arranged_curves_TM = ArrangeAntiCrossingData(MCParams{1,ii}.Detuning.E_min_MC,MCParams{1,ii}.Detuning.E_min_MCQW_TM);
            for (jj=1:length(arranged_curves_TE))
                if (~isempty(arranged_curves_TE{jj}))
                    figure(h_AntiCrossing_TE);
                    subplot(ceil(length(MCParams)/3), 3, ii); hold on; box on;
                    plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MCQW, 'b.', 'MarkerSize', 4);
                    plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MC, 'r');
                    %xlabel('E_{MC} (eV)'); ylabel('E (eV)');
                    title(['\gamma=' strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
                    %axis([1.52,1.54,1.52,1.54]);
                end
            end
            for (jj=1:length(arranged_curves_TM))
                if (~isempty(arranged_curves_TM{jj}))
                    figure(h_AntiCrossing_TM);
                    subplot(ceil(length(MCParams)/3), 3, ii); hold on; box on;
                    plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MCQW, 'b.', 'MarkerSize', 4);
                    plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MC, 'r');
                    %xlabel('E_{MC} (eV)'); ylabel('E (eV)');
                    title(['\gamma=' strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
                    %axis([1.52,1.54,1.52,1.54]);
                end
            end
            
            figure(h_AntiCrossing_TE_Waterfall);
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            PlotWaterfall(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, MCParams{1,ii}.Detuning.r_MCQW_TE, h, 'b', grid_range, y_var, 5);
            title(['\gamma=' strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_AntiCrossing_TM_Waterfall);
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            %h = subplot(1,length(MCParams), ii); hold on; box on;
            PlotWaterfall(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, MCParams{1,ii}.Detuning.r_MCQW_TM, h, 'b', grid_range, y_var, 5);
            title(['\gamma=' strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_AntiCrossing_Waterfall);  hold on;
            h = subplot(ceil(length(MCParams)/5), 5, ii); hold on; box on;
            %h = subplot(1,length(MCParams), ii); hold on; box on;
            PlotWaterfall(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, MCParams{1,ii}.Detuning.r_MCQW_TE, h, 'b', grid_range, y_var, 5);
            PlotWaterfall(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, MCParams{1,ii}.Detuning.r_MCQW_TM, h, 'r', grid_range, y_var, 5);
            title(['\gamma=' strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
            set(gca, 'YTickLabel', '');
            
            figure(h_Susceptability_TE);
            subplot(length(MCParams(1,:)),1,ii);
            [AX,H1,H2] = plotyy(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, real(MCParams{1,ii}.Detuning.n_QW_TE), MCParams{1,ii}.Detuning.E_vec/Consts.e_0, imag(MCParams{1,ii}.Detuning.n_QW_TE));
            %set(get(AX(2), 'Ylabel'), 'String', '\Im[n_{QW}]');
            set(H1,'LineStyle','-', 'Color', 'b');
            set(H2,'LineStyle','-', 'Color', 'r');
            set(AX(1), 'XLim', [min(MCParams{1,ii}.Detuning.E_vec/Consts.e_0), max(MCParams{1,ii}.Detuning.E_vec/Consts.e_0)], 'YTick', [3.5, 4.5]);
            set(AX(2), 'XLim', [min(MCParams{1,ii}.Detuning.E_vec/Consts.e_0), max(MCParams{1,ii}.Detuning.E_vec/Consts.e_0)], 'YTick', [0, 0.5]);
            set(AX(1), 'YLim', [3.5, max(real(MCParams{1,1}.Detuning.n_QW_TE))]);
            set(AX(2), 'YLim', [0, max(imag(MCParams{1,1}.Detuning.n_QW_TE))]);
            set(AX(1), 'YColor', 'b');
            set(AX(2), 'YColor', 'r');
            con_text = strrep([num2str(Params.gamma_vec(ii), '%1.0e')], 'e+0', 'x10^{');
            con_text = [con_text '}'];
            ylabel(con_text);
            h_curr = gca;
            if (ii==1)
                title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
                title(title_text);
            end
            if (ii~=length(MCParams(1,:)))
                set(AX(1), 'XTickLabel', '');
                set(AX(2), 'XTickLabel', '');
            end
            
            figure(h_Susceptability_TM);
            subplot(length(MCParams(1,:)),1,ii);
            [AX,H1,H2] = plotyy(MCParams{1,ii}.Detuning.E_vec/Consts.e_0, real(MCParams{1,ii}.Detuning.n_QW_TM), MCParams{1,ii}.Detuning.E_vec/Consts.e_0, imag(MCParams{1,ii}.Detuning.n_QW_TM));
            %set(get(AX(2), 'Ylabel'), 'String', '\Im[n_{QW}]');
            axis auto; axis tight;
            set(H1,'LineStyle','-', 'Color', 'b');
            set(H2,'LineStyle','-', 'Color', 'r');
            set(AX(1), 'XLim', [min(MCParams{1,ii}.Detuning.E_vec/Consts.e_0), max(MCParams{1,ii}.Detuning.E_vec/Consts.e_0)], 'YTick', [3.5, 4.5]);
            set(AX(2), 'XLim', [min(MCParams{1,ii}.Detuning.E_vec/Consts.e_0), max(MCParams{1,ii}.Detuning.E_vec/Consts.e_0)], 'YTick', [0, 0.5]);
            set(AX(1), 'YLim', [3.5, max(real(MCParams{1,1}.Detuning.n_QW_TM))]);
            set(AX(2), 'YLim', [0, max(imag(MCParams{1,1}.Detuning.n_QW_TM))]);
            set(AX(1), 'YColor', 'b');
            set(AX(2), 'YColor', 'r');
            con_text = strrep([num2str(Params.gamma_vec(ii), '%1.0e')], 'e+0', 'x10^{');
            con_text = [con_text '}'];
            ylabel(con_text);
            h_curr = gca;
            if (ii==1)
                title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
                title(title_text);
            end
            if (ii~=length(MCParams(1,:)))
                set(AX(1), 'XTickLabel', '');
                set(AX(2), 'XTickLabel', '');
            end
            
            n_QW_TE_Waterfall(ii,:) = MCParams{1,ii}.Detuning.n_QW_TE;
            n_QW_TM_Waterfall(ii,:) = MCParams{1,ii}.Detuning.n_QW_TM;
        end
        figure(h_Susceptability_TE);
        xlabel('E (eV)');
        figure(h_Susceptability_TM);
        xlabel('E (eV)');
        
        y_var.name = '\gamma';
        y_var.value = Params.gamma_vec;
        y_var.units = 'sec^{-1}';
        grid_range = [1.52, 1.54];
        
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, real(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'b', grid_range, y_var, 10);
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['QW Refractive Index (TE): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');
        
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, real(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'b', grid_range, y_var, 10);
        PlotWaterfall(MCParams{1,1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['QW Refractive Index (TM): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');
        
    case 3,    % single DEG concentration and single gamma value
        
        y_var.name = '\delta';
        y_var.value = nan; %  MCParams{1,1}.Detuning.delta_vec;
        y_var.units = '';
        grid_range = [1.51, 1.55];
        
        arranged_curves_TE = ArrangeAntiCrossingData(MCParams{1}.Detuning.E_min_MC,MCParams{1}.Detuning.E_min_MCQW_TE);
        arranged_curves_TM = ArrangeAntiCrossingData(MCParams{1}.Detuning.E_min_MC,MCParams{1}.Detuning.E_min_MCQW_TM);
        for (jj=1:length(arranged_curves_TE))
            if (~isempty(arranged_curves_TE{jj}))
                figure(h_AntiCrossing_TE);
                hold on; box on;
                plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MCQW, 'b.', 'MarkerSize', 8);
                plot(arranged_curves_TE{jj}.E_MC, arranged_curves_TE{jj}.E_MC, 'r');
                title(['\gamma=' strrep(num2str(Params.gamma_vec(1),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
                axis([grid_range grid_range]);
            end
        end
        for (jj=1:length(arranged_curves_TM))
            if (~isempty(arranged_curves_TM{jj}))
                figure(h_AntiCrossing_TM);
                hold on; box on;
                plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MCQW, 'b.', 'MarkerSize', 8);
                plot(arranged_curves_TM{jj}.E_MC, arranged_curves_TM{jj}.E_MC, 'r');
                %xlabel('E_{MC} (eV)'); ylabel('E (eV)');
                title(['\gamma=' strrep(num2str(Params.gamma_vec(1),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
                axis([grid_range grid_range]);
            end
        end
        
        figure(h_AntiCrossing_TE_Waterfall);
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, MCParams{1}.Detuning.r_MCQW_TE, h_AntiCrossing_TE_Waterfall, 'b', grid_range, y_var, 13);
        title(['N_{2DEG}=' strrep(num2str(Params.N_DEG,'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
        set(gca, 'YTickLabel', '');
        
        figure(h_AntiCrossing_TM_Waterfall);
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, MCParams{1}.Detuning.r_MCQW_TM, h_AntiCrossing_TM_Waterfall, 'b', grid_range, y_var, 13);
        title(['N_{2DEG}=' strrep(num2str(Params.N_DEG,'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
        set(gca, 'YTickLabel', '');
        
        figure(h_AntiCrossing_Waterfall); hold on;
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, MCParams{1}.Detuning.r_MCQW_TE, h_AntiCrossing_Waterfall, 'b', grid_range, y_var, 13);
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, MCParams{1}.Detuning.r_MCQW_TM, h_AntiCrossing_Waterfall, 'r', grid_range, y_var, 13);
        title(['N_{2DEG}=' strrep(num2str(Params.N_DEG,'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
        set(gca, 'YTickLabel', '');
        
        figure(h_Susceptability_TE);
        [AX,H1,H2] = plotyy(MCParams{1}.Detuning.E_vec/Consts.e_0, real(MCParams{1}.Detuning.n_QW_TE), MCParams{1}.Detuning.E_vec/Consts.e_0, imag(MCParams{1}.Detuning.n_QW_TE));
        set(H1,'LineStyle','-', 'Color', 'b');
        set(H2,'LineStyle','-', 'Color', 'r');
        set(AX(1), 'XLim', [min(MCParams{1}.Detuning.E_vec/Consts.e_0), max(MCParams{1}.Detuning.E_vec/Consts.e_0)], 'YTick', [3.5, 4.5]);
        set(AX(2), 'XLim', [min(MCParams{1}.Detuning.E_vec/Consts.e_0), max(MCParams{1}.Detuning.E_vec/Consts.e_0)], 'YTick', [0, 0.5]);
        set(AX(1), 'YLim', [3.5, max(real(MCParams{1}.Detuning.n_QW_TE))]);
        set(AX(2), 'YLim', [0, max(imag(MCParams{1}.Detuning.n_QW_TE))]);
        set(AX(1), 'YColor', 'b');
        set(AX(2), 'YColor', 'r');
        con_text = strrep([num2str(Params.gamma_vec(1), '%1.0e')], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
        title(title_text);
        
        figure(h_Susceptability_TM);
        [AX,H1,H2] = plotyy(MCParams{1}.Detuning.E_vec/Consts.e_0, real(MCParams{1}.Detuning.n_QW_TM), MCParams{1}.Detuning.E_vec/Consts.e_0, imag(MCParams{1}.Detuning.n_QW_TM));
        axis auto; axis tight;
        set(H1,'LineStyle','-', 'Color', 'b');
        set(H2,'LineStyle','-', 'Color', 'r');
        set(AX(1), 'XLim', [min(MCParams{1}.Detuning.E_vec/Consts.e_0), max(MCParams{1}.Detuning.E_vec/Consts.e_0)], 'YTick', [3.5, 4.5]);
        set(AX(2), 'XLim', [min(MCParams{1}.Detuning.E_vec/Consts.e_0), max(MCParams{1}.Detuning.E_vec/Consts.e_0)], 'YTick', [0, 0.5]);
        set(AX(1), 'YLim', [3.5, max(real(MCParams{1}.Detuning.n_QW_TM))]);
        set(AX(2), 'YLim', [0, max(imag(MCParams{1}.Detuning.n_QW_TM))]);
        set(AX(1), 'YColor', 'b');
        set(AX(2), 'YColor', 'r');
        con_text = strrep([num2str(Params.gamma_vec(1), '%1.0e')], 'e+0', 'x10^{');
        con_text = [con_text '}'];
        ylabel(con_text);
        h_curr = gca;
        
        title_text = ['Re[n_{QW}] (blue), Im[n_{QW}] (red): T=' num2str(Params.T) 'K'];
        title(title_text);
        
        figure(h_Susceptability_TE);
        xlabel('E (eV)');
        figure(h_Susceptability_TM);
        xlabel('E (eV)');
        
        n_QW_TE_Waterfall = MCParams{1}.Detuning.n_QW_TE;
        n_QW_TM_Waterfall = MCParams{1}.Detuning.n_QW_TM;
        
        y_var.name = 'N_{2GED}';
        y_var.value = Params.N_DEG_vec;
        y_var.units = 'cm^{-2}';
        y_var.num_type = 'exp';
        grid_range = [1.3, 1.7];
        
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, real(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'b', grid_range, y_var);
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TE_Waterfall), h_Susceptability_TE_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['Refractive Index (TE): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');
        
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, real(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'b', grid_range, y_var);
        PlotWaterfall(MCParams{1}.Detuning.E_vec/Consts.e_0, imag(n_QW_TM_Waterfall), h_Susceptability_TM_Waterfall, 'r', grid_range);
        ylabel('Re[n_{QW}] (blue), Im[n_{QW}] (red) (a.u.)');
        xlabel('E (eV)');
        title(['Refractive Index (TM): T=' num2str(Params.T) 'K']);
        set(gca, 'YTickLabel', '');
               
end

function f = SortMatrix(f)

[len_x,len_y] = size(f);
if (len_x > len_y)
    f = f.';
end

for (ii=1:len_y)
    f(:,ii) = sort(f(:,ii));
end
