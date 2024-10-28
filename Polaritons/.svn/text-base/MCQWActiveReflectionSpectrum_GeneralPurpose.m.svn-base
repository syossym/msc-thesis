clc; clear all; close all; warning off;

%% Init Simulation

global Consts;
global R G_TE G_TM z_exp;
global project_path;

project_path = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code';
cd(project_path);
run('.\Common\AddPath.m');

Constants;

%% Definitions

E_vec = [1.52:0.00005:1.54]*Consts.e_0;
lambda_vec = (2*pi./E_vec)*Consts.c*Consts.hbar; %[1300:0.1:1600].*1e-9;  % [m]
grid = 300;
E_in = [1;0];     % left-hand side incident field amplitude
n_c_l = 1; %M_GaAs.n; % refractive index of the left cladding
n_c_r = 1; %M_GaAs.n; % refractive index of the right cladding

%% Creating structure

MC_Structure = ReadStructureFile();
delta = 0.95574;

h = waitbar(0,'Building Structure...');
for (ss=1:length(MC_Structure))
    waitbar(ss/length(MC_Structure),h);
    params.x = MC_Structure{ss}.x;
    params.E = 1.525; %(E_vec(1)+0.5*(E_vec(end)-E_vec(1)))/Consts.e_0;
    mat = GetMaterial(MC_Structure{ss}.Name, params);
    n_vec_profile(ss) = mat.n;
    n_vec_calc_init(ss,:) = GetRefractiveIndex(MC_Structure{ss}.Name, E_vec./Consts.e_0);
    l_vec(ss) = delta*MC_Structure{ss}.L*1e-10;      % [m]
end
close(h);

z_vec = 0;
z_grid_calc = 0;
n_profile = n_vec_profile(1);
for (ll=1:length(l_vec))
    z_vec = [z_vec, z_vec(end)+l_vec(ll)];
    vec = linspace(z_grid_calc(end), z_grid_calc(end)+l_vec(ll), grid);
    z_grid_calc = [z_grid_calc, vec ];
    n_profile = [n_profile, ones(1,length(vec))*n_vec_profile(ll)];
end
z_grid_calc = z_grid_calc(2:end);
n_profile = n_profile(2:end);

% Load the QW simulation results
[filename, pathname, filterindex] = uigetfile({'*.mat', 'Matlab Files'; '*.*',  'All Files'});
load([pathname filename]);

%% Reflection spectrum - MC

h = waitbar(0,'Calculating MC Spectrum');
for (gg=1:length(lambda_vec))
    waitbar(gg/length(lambda_vec),h);
  
    n_vec = [n_c_l; n_vec_calc_init(:,gg)];
    %n_vec = [n_c_l, n_vec_profile];
    M_DBR = eye(2);
    M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                               (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
    M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                            (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
     
    %M_DBR = M_c_l*M_DBR;
    %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
    for (nn=1:length(l_vec))
        k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
        M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
                         0           , exp(-1i*k*l_vec(nn))];
        M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
                                    (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];
        
        M_DBR = M_p*M_i*M_DBR;
    end
    %M_DBR = M_c_r*M_DBR;
    r_MC(gg) = -M_DBR(2,1)/M_DBR(2,2);
end
close(h);

%% Reflection spectrum - MC+QW

params.x = 1;
params.E = E_vec/Consts.e_0;
QW_Mat = GetMaterial('GaAs', params);
Xi_QW = (QW_Mat.n.^2-1)/(4*pi);

r = zeros(length(con_vec), length(E_vec));
for (con_num=1:length(Xi_TE_HF_vec))
    Xi_e_Re = interp1(E_exc, smooth(real(Xi_TE_HF_vec{con_num})), E_vec, 'pchip');
    Xi_e_Im = interp1(E_exc, smooth(imag(Xi_TE_HF_vec{con_num})), E_vec, 'pchip');
    Xi_QW = Xi_QW + Xi_e_Re + 1i*Xi_e_Im;
    n_QW(con_num, :) = sqrt(1+4*pi*Xi_QW);
    n_vec_calc = n_vec_calc_init;
    
    for (ss=1:length(MC_Structure))
        if (MC_Structure{ss}.Active)
            n_vec_calc(ss,:) = n_QW(con_num, :);
        end
    end
    
    h = waitbar(0,['Calculating MCQW Spectrum, N_{DEG}=' num2str(con_vec(con_num),'%1.0e') 'cm^{-2}']);
    for (gg=1:length(lambda_vec))
        waitbar(gg/length(lambda_vec),h);
        n_vec_active = [n_c_l; n_vec_calc(:,gg)];
                
        %n_vec_active = n_vec + [0;dn_vec(:,gg)];
    
        %n_vec = [n_c_l, n_vec_profile];
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec_active(1))).*[(n_vec_active(1)+n_c_l), (n_vec_active(1)-n_c_l);...
                                   (n_vec_active(1)-n_c_l), (n_vec_active(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec_active(end)), (n_c_r-n_vec_active(end));...
                                (n_c_r-n_vec_active(end)), (n_c_r+n_vec_active(end))];

        %M_DBR = M_c_l*M_DBR;
        %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
        for (nn=1:length(l_vec))
            k = (2*pi/lambda_vec(gg))*n_vec_active(nn+1);
            M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
                             0           , exp(-1i*k*l_vec(nn))];
            M_i = (1/(2*n_vec_active(nn+1))).*[(n_vec_active(nn+1)+n_vec_active(nn)) , (n_vec_active(nn+1)-n_vec_active(nn)) ; ...
                                        (n_vec_active(nn+1)-n_vec_active(nn)) , (n_vec_active(nn+1)+n_vec_active(nn))];

            M_DBR = M_p*M_i*M_DBR;
        end
        %M_DBR = M_c_r*M_DBR;
        r_MCQW(con_num,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    close(h);
end

% Plotting
figure(2);
subplot(211);
%plot(l_profile(1:length(n_profile))./1e-9, n_profile); hold on;
plot(z_grid_calc/1e-10, n_profile); hold on;
axis([min(z_grid_calc)/1e-10, max(z_grid_calc)/1e-10, 0 , max(n_profile)+0.2]);
xlabel('z [A]'); ylabel('n');
title('(a)');
subplot(212); box on; hold on;
plot(E_vec/Consts.e_0, smooth(abs(r_MC).^2), 'r');
for (con_num=1:length(Xi_TE_HF_vec))
    plot(E_vec/Consts.e_0, smooth(abs(r_MCQW(con_num,:)).^2));
    [min_v, min_i] = min(smooth(abs(r_MCQW(con_num,:)).^2));
    text(E_vec(min_i)/Consts.e_0, min_v, [strrep(num2str(con_vec(con_num),'%1.0e'), 'e+0', 'x10^{') '}'], 'FontSize', 7);
end
title('(b)');
xlabel('E [eV]');
ylabel('|r_{MC}|^2,|r_{MCQW}|^2');

figure(3);
for (con_num=1:length(Xi_TE_HF_vec))
    [maxs_r, mins_r] = peakdet(abs(r_MCQW(con_num,:)).^2, 0.1, E_vec/Consts.e_0);
    
    subplot(length(Xi_TE_HF_vec),1,con_num);
    plot(E_vec/Consts.e_0, abs(r_MCQW(con_num,:)).^2, 'b', E_vec/Consts.e_0, abs(r_MC).^2, 'r');
    axis([min(E_vec)/Consts.e_0, max(E_vec)/Consts.e_0, 0.25, 1]);
    for (mn=1:length(mins_r(:,1)))
        text(mins_r(mn,1), mins_r(mn,2), [num2str(mins_r(mn,1)) 'eV']);
    end
    con_text = strrep([num2str(con_vec(con_num), '%1.0e')], 'e+0', 'x10^{');
    con_text = [con_text '}'];
    ylabel(con_text);
    h_curr = gca;
    if (con_num==1)
        [min_MC,min_MC_i] = min(abs(r_MC));
        title_text = ['|r_{MC}|^2 (red), |r_{MCQW}|^2 (blue): E_{MC}=' num2str(E_vec(min_MC_i)/Consts.e_0) 'eV, T=' num2str(T) 'K, \delta=' num2str(delta)];
        title(title_text);
    end
    if (con_num~=length(Xi_TE_HF_vec))
        set(h_curr, 'XTickLabel', '');
    end
end
xlabel('E [eV]');

figure(4);
for (con_num=1:length(Xi_TE_HF_vec))
    subplot(length(Xi_TE_HF_vec),1,con_num);
    [AX,H1,H2] = plotyy(E_vec/Consts.e_0, real(n_QW(con_num, :)), E_vec/Consts.e_0, imag(n_QW(con_num, :)));
    %set(get(AX(2), 'Ylabel'), 'String', '\Im[n_{QW}]');
    set(H1,'LineStyle','-', 'Color', 'b');
    set(H2,'LineStyle','-', 'Color', 'r');
    set(AX(1), 'XLim', [min(E_vec/Consts.e_0), max(E_vec/Consts.e_0)], 'YTick', [3.5, 4.5]);
    set(AX(2), 'XLim', [min(E_vec/Consts.e_0), max(E_vec/Consts.e_0)], 'YTick', [0, 0.5]);
    set(AX(1), 'YLim', [3.5, 4.5]);
    set(AX(2), 'YLim', [0, 0.5]);
    set(AX(1), 'YColor', 'b');
    set(AX(2), 'YColor', 'r');
    con_text = strrep([num2str(con_vec(con_num), '%1.0e')], 'e+0', 'x10^{');
    con_text = [con_text '}'];
    ylabel(con_text);
    h_curr = gca;
    if (con_num==1)
        title_text = ['\Re[n_{QW}] (blue), \Im[n_{QW}] (red): T=' num2str(T) 'K, \delta=' num2str(delta)];
        title(title_text);
    end
    if (con_num~=length(Xi_TE_HF_vec))
        set(AX(1), 'XTickLabel', '');
        set(AX(2), 'XTickLabel', '');
    end
end
xlabel('E [eV]');

%% Reflection spectrum - MC+QW - detuning

con_num = 7;
num_mins = 4;
N_DEG = con_vec(con_num);
delta_vec = [0.945:0.0002:0.958];

params.x = 1;
params.E = E_vec/Consts.e_0;
QW_Mat = GetMaterial('GaAs', params);
Xi_QW = (QW_Mat.n.^2-1)/(4*pi);
Xi_e_Re = interp1(E_exc, (real(Xi_TE_HF_vec{con_num})), E_vec, 'pchip');
Xi_e_Im = interp1(E_exc, (imag(Xi_TE_HF_vec{con_num})), E_vec, 'pchip');
Xi_QW = Xi_QW + Xi_e_Re + 1i*Xi_e_Im;
n_vec_calc = n_vec_calc_init;

for (ss=1:length(MC_Structure))
    if (MC_Structure{ss}.Active)
        n_vec_calc(ss,:) = sqrt(1+4*pi*Xi_QW);
    end
end

l_vec_init = l_vec/delta;
E_min_r_MCQW = zeros(num_mins, length(delta_vec));
for (dd=1:length(delta_vec))
    delta_vec(dd)
    l_vec = l_vec_init.*delta_vec(dd);
    
    %h = waitbar(0,['Calculating MC Spectrum, \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_init(:,gg)];
        %n_vec = [n_c_l, n_vec_profile];
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                                   (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                                (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];

        %M_DBR = M_c_l*M_DBR;
        %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
        for (nn=1:length(l_vec))
            k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
            M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
                             0           , exp(-1i*k*l_vec(nn))];
            M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
                                        (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];

            M_DBR = M_p*M_i*M_DBR;
        end
        %M_DBR = M_c_r*M_DBR;
        r_MC_detuning(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);
    
    [min_r_MC, min_r_MC_i] = min(abs(r_MC_detuning(dd,:)).^2);
    E_min_r_MC(dd) = E_vec(min_r_MC_i)/Consts.e_0;
    
    %h = waitbar(0,['Calculating MCQW Spectrum, \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc(:,gg)];
        %n_vec = [n_c_l, n_vec_profile];
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                                   (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                                (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];

        %M_DBR = M_c_l*M_DBR;
        %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
        for (nn=1:length(l_vec))
            k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
            M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
                             0           , exp(-1i*k*l_vec(nn))];
            M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
                                        (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];

            M_DBR = M_p*M_i*M_DBR;
        end
        %M_DBR = M_c_r*M_DBR;
        r_MCQW_detuning(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);
       
    %[maxs_r_MCQW,mins_r_MCQW] = GetExtrema(E_vec/Consts.e_0, (abs(r_MCQW_detuning(dd,:)).^2), num_mins, figure(6));
    [maxs_r_MCQW,mins_r_MCQW] = peakdet(smooth(abs(r_MCQW_detuning(dd,:)).^2), 0.01, E_vec/Consts.e_0);
    %[min, min_i, max, max_i] = extrema(abs(r_MCQW_detuning(dd,:)).^2);
    E_min_r_MCQW(1:length(mins_r_MCQW(:,1)),dd) = mins_r_MCQW(:,1);
%     for (mn=1:length(mins_r_MCQW(:,1)))
%        E_min_r_MCQW{mn}(dd) = mins_r_MCQW(mn,1); 
%     end  
end

% Plotting
figure(5);
for (dd=1:1:length(delta_vec)) 
    [maxs_r,mins_r] = peakdet(abs(r_MCQW_detuning(dd,:)).^2, 0.001, E_vec/Consts.e_0);
    subplot(length(delta_vec),1,dd);
    plot(E_vec/Consts.e_0, abs(r_MCQW_detuning(dd,:)).^2, 'b', E_vec/Consts.e_0, abs(r_MC_detuning(dd,:)).^2, 'r');
    axis([min(E_vec)/Consts.e_0, max(E_vec)/Consts.e_0, 0.25, 1]);
    ylabel(num2str(delta_vec(dd)), 'FontSize', 7);
    for (mn=1:length(mins_r(:,1)))
        text(mins_r(mn,1), mins_r(mn,2), [num2str(mins_r(mn,1)) 'eV']);
    end
    h_curr = gca;
    if (dd==1)
        title_text = ['|r_{MC}|^2 (red), |r_{MCQW}|^2 (blue): T=' num2str(T) 'K, N_{DEG}=' strrep(num2str(N_DEG, '%1.0e'), 'e+0', 'x10^{') '}cm^{-2}'];
        title(title_text);
    end
    if (dd~=length(delta_vec))
        set(h_curr, 'XTickLabel', '');
    end
end
xlabel('E [eV]');

figure(6); box on; hold on;
for (dd=1:1:length(delta_vec)) 
    plot(E_vec/Consts.e_0, abs(r_MCQW_detuning(dd,:)).^2, 'b');
    [min_E_MCQW, min_E_MCQW_i] = min(abs(r_MCQW_detuning(dd,:)).^2);
    text(E_vec(min_E_MCQW_i)/Consts.e_0, min_E_MCQW, num2str(delta_vec(dd)), 'FontSize', 7);
end
ylabel('|r_{MCQW}|.^2'); xlabel('E [eV]');

%% Anticrossing diagram fitting (coupled oscillator model)
   
params_0(1) = 1.523;
params_0(2) = 1.525;
params_0(3) = 1.7e-3;
params_0(4) = 1e-3;
params_0(5) = 1e-3;
params_0(6) = 0.6e-3;
params_0(7) = 0.7e-3;

data = sort(E_min_r_MCQW(1:3,:), 1);
%data = sort([E_min_r_MCQW{1}; E_min_r_MCQW{2}; E_min_r_MCQW{3}]);
options = optimset('TolFun',1e-1,'TolX',1e-1,'MaxIter',5000000,'MaxFunEvals',5000000, 'Display', 'iter');

params_fit = lsqcurvefit(@CoupledOscillatorModelFunction,params_0,E_min_r_MC, data);
r_MCQW_fit = CoupledOscillatorModelFunction(params_fit, E_vec/Consts.e_0);

% Plotting
h_7 = figure(7); box on; hold on; 
for (ii=1:length(data(:,1)))
    h_7_1 = plot(E_vec/Consts.e_0, r_MCQW_fit(ii,:));  
    h_7_2 = plot(E_min_r_MC, data(ii,:), '.r');
end
fit_params_text(1) = {strcat('$$E_{X_1}=',num2str(params_fit(1)), 'eV$$')};
fit_params_text(2) = {strcat('$$E_{X_2}=', num2str(params_fit(2)), 'eV$$')};
fit_params_text(3) = {strcat('$$\hbar\Omega_{X_1}=', num2str(params_fit(3)/1e-3), 'meV$$')};
fit_params_text(4) = {strcat('$$\hbar\Omega_{X_2}=', num2str(params_fit(4)/1e-3), 'meV$$')};
fit_params_text(5) = {strcat('$$\gamma_{MC}=', num2str(params_fit(5)), 'eV$$')};
fit_params_text(6) = {strcat('$$\gamma_{X_1}=', num2str(params_fit(6)), 'eV$$')};
fit_params_text(7) = {strcat('$$\gamma_{X_2}=', num2str(params_fit(7)), 'eV$$')};
text(max(E_vec)/Consts.e_0,max(E_vec)/Consts.e_0,fit_params_text,'HorizontalAlignment','left','interpreter','latex' );
xlabel('$$E_{MC} [eV]$$','interpreter','latex'); ylabel('$$E [eV]$$','interpreter','latex');
title(['Level anticrossing diagram: $$T=' num2str(T) 'K$$, $$N_{2DEG}=' strrep([num2str(con_vec(con_num), '%1.0e')], 'e+0', '\times 10^{') '}cm^{-2}$$'],'interpreter','latex');
axis tight; set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
