function HartreeFock_Test2()

global Consts;
global R G_TE G_TM z_exp;

%% Part 1 - Subband energies
% close all; clc; %clear all;
% warning off;
% 
% %Init project
% global project_path;
% project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
% cd(project_path);
% run('.\Common\AddPath.m');
% 
% Constants;

% well_width = 50   % [A]
% x = 0.2;
% Res1 = SingleQW_kp(well_width,x);

L_z = 200e-10;

% Structure = { 'GaAlAs' , 300, 0.1 ;
%     'GaAs', L_z*1e10, 0;
%     'GaAlAs' , 300, 0.1 };
% 
% Res2 = GeneralStructure_TransferMatrix(Structure);
% 
% R = Res2;
% 
% close all;

%% Part 2 - Gain Spectra

% Definitions
con_vec = [2e11,6e11,2e12,3e12,4e12,6e12,8e12];    % [cm^-2]
T = 300;
gamma = 1e13;                                      % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (R.Materials{2}.E_g+0.5)*Consts.e_0;
Mb_s = Consts.e_0*R.Materials{2}.E_p*Consts.m_0/6;
C_g = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
init_E = 0.01*Consts.e_0;
cb_len = 1; vb_len = 3;

% Defining simulation grids
% z_exp = exp(-abs(repmat(R.z_grid,1,length(R.z_grid)) -...
%                  repmat(R.z_grid',length(R.z_grid),1)));
sim_index = 2:350;
k_vec = R.k_t_mat(sim_index);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.1:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) -...
                 repmat(R.z_grid',length(R.z_grid),1));

% for (tt = 1:length(theta))
%     q_mat(:,:,tt) = sqrt(k_mat.^2 + k_mat.'.^2 - 2.*k_mat.*k_mat.'.*cos(theta(tt)));
% end

for (con_num = 1:length(con_vec))
    
    % Carrier concentrations
    N = con_vec(con_num);
    P = N;
    
    % Quasi Fermi levels
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
    E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);
    
    for (ee = 1:length(E_exc))
        E_exc(ee)
        w = (E_exc(ee)./Consts.hbar);
        sum_TE = 0;
        sum_TM = 0;
        
        for (cb_num = 1:cb_len)
            for (vb_num = 1:vb_len)
                disp(['cb=' num2str(cb_num) ', vb_num=' num2str(vb_num)]);

                % Carrier distributions -----------------------------------
                E_e = R.E_k_c{cb_num}(sim_index);
                E_h = R.E_k_v{vb_num}(sim_index);
                f_h = FermiDirac(E_fv, E_h, T);
                f_e = FermiDirac(E_fc, E_e, T);
                omega = f_e + f_h - 1;
                
                % Simulation definitions ----------------------------------
                mu_TE = sqrt(squeeze(R.rtsTE(vb_num,cb_num,sim_index))).'.*Mb_s;
                mu_TM = sqrt(squeeze(R.rtsTM(vb_num,cb_num,sim_index))).'.*Mb_s;
                Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
                Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
                Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
                Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
                Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
                Xi_ee = Xi_e.*Xi_e.';
                Xi_hh = Xi_h.'.*Xi_h;
                Xi_eh = Xi_e.*Xi_h;
                
                % Theta matrix calculation --------------------------------
                disp('-- Calculating Theta matrix:');
                C_Theta = Consts.e_0^2./(8*pi^2*R.Materials{2}.eps_r*Consts.eps_0);
                Theta = zeros(k_len);
                for (ii = 1:k_len)
                    for (jj = 1:k_len)
                        disp(['i=' num2str(ii) ', j=' num2str(jj)]);
                        q = sqrt(k_vec(ii).^2 + k_vec(jj).^2 - 2.* k_vec(ii).*k_vec(jj).*cos(theta));
                        q(q==0) = 1e-5;
                        Gq_eh = zeros(1,length(theta));
                        eps_q = zeros(1,length(theta));
                        for (tt = 1:length(theta))
                            % 2D form factor - hole-electron
                            Gq_eh(tt) = G_q(Xi_eh,q(tt),R.z_grid);
                            % Dielectric function
                            eps_q(tt) = Eps_q(k_vec,q(tt),f_e,f_h,E_e,E_h,theta,R.Materials{2}.eps_r);
                        end
                        
                        Theta(ii,jj) = C_Theta.*trapz(theta, Gq_eh./(q.*eps_q));
                        
                        figure(1); hold on;
                        subplot(211); plot(q,Gq_eh);
                        ylabel('G_q');
                        subplot(212); plot(q,eps_q);
                        xlabel('q'); ylabel('\epsilon_q');
                        drawnow;
                    end
                end
                
                % Xi matrix calculation -----------------------------------
                disp('-- Calculating D_E_CH:');
                q_vec = k_vec; q_vec(1) = 1e-5;
                D_E_CH = E_CH(Xi_hh,R.z_grid,k_vec,q_vec,f_e,f_h,E_e,E_h,theta,R.Materials{2}.eps_r);
                
                % Screened-exchange shift energy
                disp('-- Calculating D_E_SX:');
                D_E_SX = E_SX(Xi_ee,Xi_hh,R.z_grid,k_vec,q_vec,f_e,f_h,E_e,E_h,theta,R.Materials{2}.eps_r);
                
                % Bandgap renormalization
                w_tag = (E_h + E_e + D_E_CH + D_E_SX)./Consts.hbar;
                
                Xi_0_TE = (-1i/Consts.hbar).*((omega.*mu_TE)./(1i*(w_tag-w)+gamma));
                Xi_0_TM = (-1i/Consts.hbar).*((omega.*mu_TM)./(1i*(w_tag-w)+gamma));
                
                % Coulomb enhancement factor calculation ------------------
                M_TE = eye(length(R.k_t_mat));
                M_TE = M_TE - repmat((dk./mu_TE).*Xi_0_TE, k_len, 1).*Theta;
                M_TM = eye(length(R.k_t_mat));
                M_TM = M_TM - repmat((dk./mu_TM).*Xi_0_TM, k_len, 1).*Theta;
                Q_k_TE = M_TE\ones(k_len,1);
                Q_k_TM = M_TM\ones(k_len,1);
                
                % Calculating gain sum elements ---------------------------
                I_TE = 2*trapz(R.k_t_mat, k_vec.*mu_TE^2.*omega.*(1./((w_tag-w).^2+gamma^2)).*Q_k_TE);
                I_TM = 2*trapz(R.k_t_mat, k_vec.*mu_TM^2.*omega.*(1./((w_tag-w).^2+gamma^2)).*Q_k_TM);
                sum_TE = sum_TE + I_TE;
                sum_TM = sum_TM + I_TM;
            end
        end
        
        G_TE(ee) = imag(C.*1i*w.*sum_TE*L_z);
        G_TM(ee) = imag(C.*1i*w.*sum_TM*L_z);
        
    end
    
    index_TE = find(abs(G_TE) == max(abs(G_TE(1:round(length(G_TE)/3)))));
    index_TM = find(abs(G_TM) == max(abs(G_TM(1:round(length(G_TM)/3)))));
    
    figure(1);
    subplot(211);
    plot(E_exc/Consts.e_0, G_TE/100); hold on;
    plot(Res2.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
    ylabel('G [cm^{-1}]');
    text(E_exc(index_TE)/Consts.e_0,G_TE(index_TE)/100, [num2str(N/1e10), 'x10^{10} cm^{-1}'], 'FontSize', 8);
    subplot(212);
    plot(E_exc/Consts.e_0, G_TM/100); hold on;
    plot(Res2.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
    text(E_exc(index_TM)/Consts.e_0,G_TM(index_TM)/100, [num2str(N/1e10), 'x10^{10} cm^{-1}'], 'FontSize', 8);
    xlabel('E [eV]'); ylabel('G [cm^{-1}]');
end

return;

%
% This function calculates the 2D form factor
%
function G = G_q(Xi12,q,z_grid)

global z_exp;
z_exp_q = exp(z_exp.*q);
integrand_G = Xi12.*z_exp_q;
G = trapz(z_grid, trapz(z_grid,integrand_G,1), 2);

%
% This function calculates the screened dielectric function
%
function eps = Eps_q(k,q,f_e,f_h,E_e,E_h,theta,eps_r)

global Consts;
integrand_e = zeros(size(k));
integrand_h = integrand_e;
for(kk = 1:length(k))
    k_minus_q = sqrt(k(kk).^2+q.^2-2.*k(kk).*q.*cos(theta));
    f_e_k_minus_q = interp1(k,f_e,k_minus_q);
    f_h_k_minus_q = interp1(k,f_h,k_minus_q);
    E_e_k_minus_q = interp1(k,E_e,k_minus_q);
    E_h_k_minus_q = interp1(k,E_h,k_minus_q);
    
    integrand_e(kk) = trapz(theta, (f_e_k_minus_q - f_e(kk).*ones(size(theta)))./(E_e_k_minus_q - E_e(kk).*ones(size(theta))));
    integrand_h(kk) = trapz(theta, (f_h_k_minus_q - f_h(kk).*ones(size(theta)))./(E_h_k_minus_q - E_h(kk).*ones(size(theta))));
    
    integrand_e(isnan(integrand_e)) = 0; integrand_h(isnan(integrand_h)) = 0;
    integrand_e(integrand_e == -inf) = -1e19; integrand_h(integrand_h == -inf) = -1e19;
    
%     figure(2); 
%     plot(theta, (f_e_k_minus_q - f_e(kk))./(E_e_k_minus_q - E_e(kk)), 'b',...
%          theta, (f_h_k_minus_q - f_h(kk))./(E_h_k_minus_q - E_h(kk)), 'r');
%     xlabel('k'); ylabel('(f^e_{|k-q|} - f^e_{|k|})/(E^e_{|k-q|} - E^e_{|k|})');
%     drawnow;
end
eps = 1 - (Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0))*(trapz(k, k.*integrand_e) + trapz(k, k.*integrand_h));

%
% This function calculates the screened-exchange shift energy
%
function D_E_SX = E_SX(Xi_ee,Xi_hh,z_grid,k,q,f_e,f_h,E_e,E_h,theta,eps_r)

global Consts;
C_E_SX = Consts.e_0^2./(8*pi^2*R.Materials{2}.eps_r*Consts.eps_0);
D_E_SX = zeros(1,length(k)); int_E_SX = zeros(1,length(k));
for (kk_tag = 1:length(k))
    for (kk = 1:length(k))
        k_minus_k_tag = sqrt(k(kk_tag).^2+k(kk).^2-2.*k(kk_tag).*k(kk).*cos(theta));
        int_E_SX(kk) = trapz(theta, (G_q(Xi_ee,k_minus_k_tag,z_grid).*f_e(kk)+...
                                     G_q(Xi_hh,k_minus_k_tag,z_grid).*f_h(kk)).*...
                                     (1./(k_minus_k_tag.*Eps_q(k,k_minus_k_tag,f_e,f_h,E_e,E_h,theta,eps_r))));
    end
    D_E_SX(kk_tag) = -C_E_SX.*trapz(k_vec, k_vec.*int_E_SX);
end

%
% This function calculates the Coulomb-hole self-energy
%
function D_E_CH = E_CH(Xi_hh,z_grid,k,q,f_e,f_h,E_e,E_h,theta,eps_r)

global Consts;
int_E_CH = zeros(1,length(q));
for (qq = 1:length(q))
    int_E_CH(qq) =  G_q(Xi_hh,q(qq),z_grid).*((1./Eps_q(k,q(qq),f_e,f_h,E_e,E_h,theta,eps_r))-1);
end
D_E_CH = Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)*trapz(q,int_E_CH);
