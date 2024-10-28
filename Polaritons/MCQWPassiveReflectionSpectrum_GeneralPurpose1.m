clc; clear all; close all; warning off;

%% Init Simulation
global Consts;
global R G_TE G_TM z_exp;
global project_path;

project_path = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code';
cd(project_path);
run('.\Common\AddPath.m');

% Create the physical constants structure
Constants;
params.x = 0.1;
params.T = 2;
params.E = 0;
M_GaAs = GetMaterial('GaAs', params);
M_AlAs = GetMaterial('AlAs', params);
M_GaAlAs = GetMaterial('GaAlAs', params);

%% Definitions

E_vec = [1.4:0.0001:1.65]*Consts.e_0;
lambda_vec = (2*pi./E_vec)*Consts.c*Consts.hbar; %[1300:0.1:1600].*1e-9;  % [m]
grid = 300;
E_in = [1;0];     % left-hand side incident field amplitude
n_c_l = 1; %M_GaAs.n; % refractive index of the left cladding
n_c_r = 1; %M_GaAs.n; % refractive index of the right cladding

%% Creating structure

Structure = ReadStructureFile();
delta = 0.96; %0.95574;

h = waitbar(0,'Building Structure...');
for (ss=1:length(Structure))
    waitbar(ss/length(Structure),h);
    params.x = Structure{ss}.x;
    params.E = 1.525; %(E_vec(1)+0.5*(E_vec(end)-E_vec(1)))/Consts.e_0;
    mat = GetMaterial(Structure{ss}.Name, params);
    n_vec_profile(ss) = mat.n;
    n_vec_calc(ss,:) = GetRefractiveIndex(Structure{ss}.Name, E_vec./Consts.e_0);
    l_vec(ss) = delta*Structure{ss}.L*1e-10;      % [m]
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

%% Simulation

% Reflection spectrum
h = waitbar(0,'Calculating Spectrum...');
I_mat = zeros(length(lambda_vec), length(l_vec));
for (gg=1:length(lambda_vec))
    waitbar(gg/length(lambda_vec),h);
  
    n_vec = [n_c_l; n_vec_calc(:,gg)];
    %n_vec = [n_c_l, n_vec_profile];
    M_DBR = eye(2);
    M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                               (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
    M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                            (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
                        
    E_mat = E_in;
           
    %M_DBR = M_c_l*M_DBR;
    %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
    for (nn=1:length(l_vec))
        k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
        M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
                         0           , exp(-1i*k*l_vec(nn))];
        M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
                                    (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];
        
        M_DBR = M_p*M_i*M_DBR;
        E_mat = M_DBR*E_mat;
        I_mat(gg,nn) = abs(E_mat(1,1)+E_mat(2,1))^2;
    end
    %M_DBR = M_c_r*M_DBR;
    r(gg) = -M_DBR(2,1)/M_DBR(2,2);
end
close(h);

% Field amplitude distribution
% M_DBR = eye(2);
% M_DBR = M_c_l*M_DBR;
% E_mat = M_DBR*E_in;
% I_mat = abs(E_mat(1,1)+E_mat(2,1))^2;
% for (nn=1:length(n_profile))
%     k = (2*pi/lambda_c)*n_profile(nn);
%     if (nn==1)
%         z = z_grid_calc(nn);
%         M =  [exp(1i*k*z),         0 ;...
%             0             , exp(-1i*k*z)];
%     elseif (n_profile(nn) == n_profile(nn-1))
%         z = z_grid_calc(nn)-z_grid_calc(nn-1);
%         M =  [exp(1i*k*z),         0 ;...
%             0             , exp(-1i*k*z)];
%     elseif (n_profile(nn) ~= n_profile(nn-1))
%         M = (1/(2*n_profile(nn))).*[(n_profile(nn)+n_profile(nn-1)) , (n_profile(nn)-n_profile(nn-1)) ; ...
%             (n_profile(nn)-n_profile(nn-1)) , (n_profile(nn)+n_profile(nn-1))];
%     end
%     M_DBR = M*M_DBR;
%     E_mat = [E_mat, M_DBR*E_in];
%     I_mat = [I_mat, abs(E_mat(1,nn)+E_mat(2,nn))^2];
% end
% M_DBR = M_c_r*M_DBR;
% E_mat = [E_mat, M_DBR*E_in];
% I_mat = [I_mat, abs(E_mat(1,end)+E_mat(2,end))^2];
% I_mat = I_mat(2:end-1);

%% Plotting

% figure(1);
% subplot(311);
% %plot(l_profile(1:length(n_profile))./1e-9, n_profile); hold on;
% plot(z_grid_calc/1e-10, n_profile); hold on;
% %vec = 0:max(n_profile)/length(n_profile):max(n_profile)-max(n_profile)/length(n_profile);
% %plot(zeros(1,length(vec)), vec, ':g',...
% %    ones(1,length(vec)).*(max(l_profile)-grid*dl)/1e-9, vec, ':g'); hold on;
% %axis([min(l_profile)./1e-9, max(l_profile)./1e-9, 0 , max(n_profile)+0.2]);
% axis([min(z_grid_calc)/1e-10, max(z_grid_calc)/1e-10, 0 , max(n_profile)+0.2]);
% xlabel('z [A]'); ylabel('n');
% title('(a)');
% subplot(312);
% %plot(lambda_vec./lambda_c, abs(r)); hold on;
% [AX,H1,H2] = plotyy(lambda_vec./lambda_c, abs(r), lambda_vec./lambda_c, angle(r)/pi);
% set(get(AX(2), 'Ylabel'), 'String', '\alpha_r/\pi');
% set(H1,'LineStyle','-', 'Color', 'b');
% set(H2,'LineStyle',':', 'Color', 'r');
% set(AX(1), 'XLim', [min(lambda_vec/lambda_c), max(lambda_vec/lambda_c)]);
% set(AX(2), 'XLim', [min(lambda_vec/lambda_c), max(lambda_vec/lambda_c)]);
% set(AX(1), 'YColor', 'k');
% set(AX(2), 'YColor', 'k');
% %plot(lambda_c./1e-9*ones(1,length(r)), abs(r), ':r');
% %text(lambda_c./1e-9, max(abs(r))-0.1, ['\lambda_c=' num2str(lambda_c/1e-9) ' nm']);
% %axis([min(lambda_vec/lambda_c), max(lambda_vec/lambda_c), 0, 1]);
% xlabel('\lambda/\lambda_c');
% ylabel('|r_{MC}|');
% title('(b)');
% subplot(313);
% plot(phi_vec/pi, abs(r)); hold on;
% text(0.5, max(abs(r))-0.1, '\phi=\pi/2');
% axis([min(phi_vec/pi), max(phi_vec/pi), min(abs(r)), 1]);
% xlabel('\phi/\pi');
% ylabel('|r_{MC}|');
% title('(c)');

figure(2);
subplot(211);
%plot(l_profile(1:length(n_profile))./1e-9, n_profile); hold on;
plot(z_grid_calc/1e-10, n_profile); hold on;
axis([min(z_grid_calc)/1e-10, max(z_grid_calc)/1e-10, 0 , max(n_profile)+0.2]);
xlabel('z [A]'); ylabel('n');
title('(a)');
subplot(212);
plot(E_vec/Consts.e_0, smooth(abs(r).^2));
xlabel('E [eV]');
ylabel('|r_{MC}|^2');

% figure(3);
% [AX,H1,H2] = plotyy(z_grid_calc/1e-10, I_mat, z_grid_calc/1e-10, n_profile);
% set(get(AX(2), 'Ylabel'), 'String', 'n');
% set(H1,'LineStyle','-', 'Color', 'r');
% set(H2,'LineStyle','-', 'Color', 'b');
% set(AX(1), 'XLim', [min(z_grid_calc/1e-10), max(z_grid_calc/1e-10)]);
% set(AX(2), 'XLim', [min(z_grid_calc/1e-10), max(z_grid_calc/1e-10)]);
% set(AX(2), 'YLim', [2, 3.5]);
% set(AX(1), 'YColor', 'k');
% set(AX(2), 'YColor', 'k');
% ylabel('|E^{+}+E^{-}|^2');
% xlabel('z [A]');