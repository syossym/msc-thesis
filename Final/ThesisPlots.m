close all;
% warning off;
%
% % Set project directories
% global project_path;
% project_path = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code';
% %project_path = 'z:\Yossi Michaeli\Documents\Thesis\Code';
% cd(project_path);
% run('.\Common\AddPath.m');
%
% % Create the physical constants structure
% Constants;

%% ==========

if (iscell(Bands))
    Temp = Bands{1};
    clear Bands;
    Bands = Temp;
end
if (iscell(QWParams))
    Temp = QWParams{1};
    clear QWParams;
    QWParams = Temp;
end

x_scale = [1.52, 1.54];
y_scale = [0, 0.01];
smooth_param = 3;

%% =========

figure(1);
subplot(211); box on; hold on;
[AX,H1,H2] = plotyy(Params.E_exc(1:length(QWParams.Xi.TE.FCT))/Consts.e_0, smooth(real(QWParams.Xi.TE.FCT),smooth_param),...
    Params.E_exc(1:length(QWParams.Xi.TE.FCT))/Consts.e_0, smooth(imag(QWParams.Xi.TE.FCT),smooth_param));
plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [-5,5], '--k');
ylabel('\Re[\chi^{TE}_{FCT}]'); title('(a)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YColor', 'k');
% set(AX(1), 'YLim', [0, 0.04]);
% set(AX(2), 'YLim', [0, 0.025]);
set(get(AX(2), 'Ylabel'), 'String', '\Im[\chi^{TE}_{FCT}]');

subplot(212); box on; hold on;
[AX,H1,H2] = plotyy(Params.E_exc(1:length(QWParams.Xi.TM.FCT))/Consts.e_0, smooth(real(QWParams.Xi.TM.FCT),smooth_param),...
    Params.E_exc(1:length(QWParams.Xi.TM.FCT))/Consts.e_0, smooth(imag(QWParams.Xi.TM.FCT),smooth_param));
plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [-5,5], '--k');
xlabel('E (eV)'); ylabel('\Re[\chi^{TM}_{FCT}]'); title('(b)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YColor', 'k', 'YTickLabelMode', 'auto');
% set(AX(1), 'YLim', [0, 0.04]);
% set(AX(2), 'YLim', [0, 0.025]);
set(get(AX(2), 'Ylabel'), 'String', '\Im[\chi^{TM}_{FCT}]');

%% =========

figure(2);
subplot(211); hold on; box on;
[AX,H1,H2] = plotyy(Params.E_exc(1:length(QWParams.Xi.TE.HF))/Consts.e_0, smooth(real(QWParams.Xi.TE.HF),smooth_param),...
    Params.E_exc(1:length(QWParams.Xi.TE.HF))/Consts.e_0, smooth(imag(QWParams.Xi.TE.HF),smooth_param));
plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [-5,5], '--k');
ylabel('\Re[\chi^{TE}_{HF}]'); title('(a)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YColor', 'k');
set(get(AX(2), 'Ylabel'), 'String', '\Im[\chi^{TE}_{HF}]');

subplot(212);
[AX,H1,H2] = plotyy(Params.E_exc(1:length(QWParams.Xi.TM.HF))/Consts.e_0, smooth(real(QWParams.Xi.TM.HF),smooth_param),...
    Params.E_exc(1:length(QWParams.Xi.TM.HF))/Consts.e_0, smooth(imag(QWParams.Xi.TM.HF),smooth_param));
hold on; box on;
%plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [-5,5], '--k');
xlabel('E (eV)'); ylabel('\Re[\chi^{TM}_{HF}]'); title('(b)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YColor', 'k', 'YTickLabelMode', 'auto');
%set(AX(1), 'YLim', [-0.4, 0.4]);
%set(AX(2), 'YLim', [0, 0.8]);
set(get(AX(2), 'Ylabel'), 'String', '\Im[\chi^{TM}_{HF}]');

%% ==========

figure(3);
subplot(211); box on; hold on;
plot(Params.E_exc(1:length(QWParams.alpha.TE.FCT))/Consts.e_0, smooth(QWParams.alpha.TE.FCT,smooth_param), 'b',...
    Params.E_exc(1:length(QWParams.alpha.TE.HF))/Consts.e_0, smooth(QWParams.alpha.TE.HF,smooth_param), 'r');
plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [0,0.5e5], '--k');
ylabel('\alpha^{TE} (cm^{-1})'); title('(a)'); legend('FCT', 'HF');
set(gca, 'XLim', x_scale)

subplot(212); box on; hold on;
plot(Params.E_exc(1:length(QWParams.alpha.TM.FCT))/Consts.e_0, smooth(QWParams.alpha.TM.FCT,smooth_param), 'b',...
    Params.E_exc(1:length(QWParams.alpha.TM.HF))/Consts.e_0, smooth(QWParams.alpha.TM.HF,smooth_param), 'r');
plot(Bands.Cond.E_k(1)./Consts.e_0*ones(1,2), [0,0.5e5], '--k');
xlabel('E (eV)'); ylabel('\alpha^{TM} (cm^{-1})'); title('(b)');
set(gca, 'XLim', x_scale);

%% ==========

figure(4);
box on; hold on;
temp_FCT = QWParams.Rsp.FCT;
temp_FCT(temp_FCT<0) = 0;
temp_HF = QWParams.Rsp.HF;
temp_HF(temp_HF<0) = 0;
plot(Params.E_exc(1:length(QWParams.Rsp.FCT))/Consts.e_0, (smooth(temp_FCT,smooth_param)), 'b',...
    Params.E_exc(1:length(QWParams.Rsp.HF))/Consts.e_0, (smooth(temp_HF,smooth_param)), 'r');
%plot(Bands.Cond.E_k(1)./Consts.e_0, [0:100:10e3], ':k');
xlabel('E (eV)'); ylabel('r_{sp} (s^{-1}m^{-3}J^{-1})'); legend('FCT', 'HF');
set(gca, 'XLim', x_scale);

%% ==========

E_vec = min(Params.E_exc):1e-5*Consts.e_0:max(Params.E_exc);
temp_params.x = 1;
temp_params.E = E_vec/Consts.e_0;
temp_params.T = Params.T;
QW_Mat = GetMaterial('GaAs', temp_params);
Xi_QW = (QW_Mat.n.^2-1)/(4*pi);
Xi_e_Re_TE_FCT = interp1(Params.E_exc, (real(QWParams.Xi.TE.FCT)), E_vec, 'pchip');
Xi_e_Im_TE_FCT = interp1(Params.E_exc, (imag(QWParams.Xi.TE.FCT)), E_vec, 'pchip');
Xi_e_Re_TM_FCT = interp1(Params.E_exc, (real(QWParams.Xi.TM.FCT)), E_vec, 'pchip');
Xi_e_Im_TM_FCT = interp1(Params.E_exc, (imag(QWParams.Xi.TM.FCT)), E_vec, 'pchip');
Xi_QW_TE_FCT = Xi_QW + Xi_e_Re_TE_FCT + 1i*Xi_e_Im_TE_FCT;
Xi_QW_TM_FCT = Xi_QW + Xi_e_Re_TM_FCT + 1i*Xi_e_Im_TM_FCT;
n_QW_TE_FCT = sqrt(1+4*pi*Xi_QW_TE_FCT);
n_QW_TM_FCT = sqrt(1+4*pi*Xi_QW_TM_FCT);
Xi_e_Re_TE_HF = interp1(Params.E_exc, (real(QWParams.Xi.TE.HF)), E_vec, 'pchip');
Xi_e_Im_TE_HF = interp1(Params.E_exc, (imag(QWParams.Xi.TE.HF)), E_vec, 'pchip');
Xi_e_Re_TM_HF = interp1(Params.E_exc, (real(QWParams.Xi.TM.HF)), E_vec, 'pchip');
Xi_e_Im_TM_HF = interp1(Params.E_exc, (imag(QWParams.Xi.TM.HF)), E_vec, 'pchip');
Xi_QW_TE_HF = Xi_QW + Xi_e_Re_TE_HF + 1i*Xi_e_Im_TE_HF;
Xi_QW_TM_HF = Xi_QW + Xi_e_Re_TM_HF + 1i*Xi_e_Im_TM_HF;
n_QW_TE_HF = sqrt(1+4*pi*Xi_QW_TE_HF);
n_QW_TM_HF = sqrt(1+4*pi*Xi_QW_TM_HF);

figure(5);
subplot(211); box on; hold on;
[AX,H1,H2] = plotyy(E_vec(1:length(n_QW_TE_HF))/Consts.e_0, smooth(real(n_QW_TE_HF),smooth_param),...
    E_vec(1:length(n_QW_TE_HF))/Consts.e_0, smooth(imag(n_QW_TE_HF),smooth_param));
plot(Bands.Cond.E_k(1)./Consts.e_0, [0:0.0001:0.05], ':k');
ylabel('\Re[n^{TE}_{HF}]'); title('(a)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YLim', [2.5, 4.5], 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YLim', [0, 1.5], 'YColor', 'k');
set(get(AX(2), 'Ylabel'), 'String', '\Im[n^{TE}_{HF}]');

subplot(212); box on; hold on;
[AX,H1,H2] = plotyy(E_vec(1:length(n_QW_TM_HF))/Consts.e_0, smooth(real(n_QW_TM_HF),smooth_param),...
    E_vec(1:length(n_QW_TM_HF))/Consts.e_0, smooth(imag(n_QW_TM_HF),smooth_param));
plot(Bands.Cond.E_k(1)./Consts.e_0, [0:0.0001:0.05], ':k');
xlabel('E (eV)'); ylabel('\Re[n^{TM}_{HF}]'); title('(b)');
set(H1,'LineStyle','-', 'Color', 'b');
set(H2,'LineStyle','-', 'Color', 'r');
set(AX(1), 'XLim', x_scale, 'YLim', [2.5, 4.5], 'YColor', 'k');
set(AX(2), 'XLim', x_scale, 'YColor', 'k', 'YLim', [0, 1.5], 'YTickLabelMode', 'auto');
set(get(AX(2), 'Ylabel'), 'String', '\Im[n^{TM}_{HF}]');

figure(6);
subplot(211);
plot(E_vec(1:length(n_QW_TE_HF))/Consts.e_0, smooth(real(n_QW_TE_HF),smooth_param));
xlabel('E (eV)'); ylabel('\Re[n^{TE}_{HF}]');
set(gca, 'XLim', x_scale);
subplot(212);
plot(E_vec(1:length(n_QW_TE_HF))/Consts.e_0, smooth(imag(n_QW_TE_HF),smooth_param));
xlabel('E (eV)'); ylabel('\Im[n^{TE}_{HF}]');
set(gca, 'XLim', x_scale);

figure(7);
subplot(211);
plot(E_vec(1:length(n_QW_TM_HF))/Consts.e_0, smooth(real(n_QW_TM_HF),smooth_param));
xlabel('E (eV)'); ylabel('\Re[n^{TM}_{HF}]');
set(gca, 'XLim', x_scale);
subplot(212);
plot(E_vec(1:length(n_QW_TM_HF))/Consts.e_0, smooth(imag(n_QW_TM_HF),smooth_param));
xlabel('E (eV)'); ylabel('\Im[n^{TM}_{HF}]');
set(gca, 'XLim', x_scale);

%% ===========================================

% father_dirname = uigetdir(project_path, 'Pick Data Directory');
% filenames = dir(father_dirname);
% names_array = [];
% names = [];
%
% for (ii=1:length(filenames))
%     if (strfind(filenames(ii).name, 'Temp_'))
%         temp1 = strsplit('_', filenames(ii).name);
%         temp2 = strcat(temp1(4), '.', temp1(5));
%         names_array = [names_array; str2double(temp2)];
%         names = [names; filenames(ii).name];
%     end
% end
% [Y,I] = sort(names_array, 'ascend');
% filenames = names(I,:);
%
% for (ii=1:length(names))

%end

%% ===========================================

cl = ['b  ';'r  ';'g  ';'c  ';'m  ';'y  ';'b  ';'r  ';'g  ';'c  ';'m  ';'y  '];

figure(8);
for (ii=1:length(QWParams))
    subplot(211); box on; hold on;
    plot(Params.E_exc./Consts.e_0, smooth(QWParams{ii}.alpha.TE.HF, smooth_param), 'Color', cl(ii));
    ylabel('\alpha_{TE}^{HF}');
    set(gca, 'XLim', x_scale);
    subplot(212); box on; hold on;
    plot(Params.E_exc./Consts.e_0, smooth(QWParams{ii}.alpha.TM.HF, smooth_param), 'Color', cl(ii));
    ylabel('\alpha_{TM}^{HF}');
    xlabel('E (eV)');
    set(gca, 'XLim', x_scale);
end

figure(9);
for (ii=1:length(QWParams))
    box on; hold on;
    plot(Params.E_exc./Consts.e_0, smooth(abs(QWParams{ii}.Rsp.HF)), 'Color', cl(ii));
    ylabel('r^{HF}_{sp}');
    xlabel('E (eV)');
    set(gca, 'XLim', x_scale);
end


%% =======================

delta_index = 15;
reflection_delta_TE = [];
reflection_delta_TM = [];
for (ii=1:length(MCParams_New))
    reflection_delta_TE = [reflection_delta_TE; MCParams_New{ii}.Detuning.HF.r_MCQW_TE(delta_index,:)];
    reflection_delta_TM = [reflection_delta_TM; MCParams_New{ii}.Detuning.HF.r_MCQW_TM(delta_index,:)];
end

grid_range = [1.52, 1.535];
y_var.name = 'N_{2GED}';
y_var.value = Params.N_DEG_vec;
y_var.units = 'cm^{-2}';
y_var.num_type = 'exp';

h_TE = figure(1);
PlotWaterfall(MCParams_New{1,1}.Detuning.HF.E_vec/Consts.e_0, abs(reflection_delta_TE), h_TE, 'b', grid_range, y_var);
title('(a)'); xlabel('E (eV)'); ylabel('Reflection (a.u.)');

h_TM = figure(2);
PlotWaterfall(MCParams_New{1,1}.Detuning.HF.E_vec/Consts.e_0, abs(reflection_delta_TM), h_TM, 'b', grid_range, y_var);
title('(b)'); xlabel('E (eV)'); ylabel('Reflection (a.u.)');

