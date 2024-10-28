global Consts;

%% Part 1 - Subband energies
% close all; clc; clear all;
% warning off;
% 
% %Init project
% global project_path;
% project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
% cd(project_path);
% run('.\Common\AddPath.m');
% 
% Constants;
% 
% % well_width = 50   % [A]
% % x = 0.2;
% % Res1 = SingleQW_kp(well_width,x);
% 
% L_z = 200e-10;
% 
% Structure = { 'GaAlAs' , 300, 0.1 ;
%               'GaAs', L_z*1e10, 0;
%               'GaAlAs' , 300, 0.1 };
% 
% Res2 = GeneralStructure_TransferMatrix(Structure);
% 
% R = Res2;
% 
% close all;

%% Part 2 - Gain Spectra

clear Gq_h Gq_l Vq_h Vq_l;
Con_vec = [2e11,6e11,2e12,3e12,4e12,6e12,8e12];    % [cm^-2]
T = 20;
gamma = 1e13;                                      % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (R.Materials{2}.E_g+0.5)*Consts.e_0;
Mb_s = Consts.e_0*Res2.Materials{2}.E_p*Consts.m_0/6;
C_g = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
init_E = 0.01*Consts.e_0;
cb_len = 1; vb_len = 3;

for (con_num = 1:length(Con_vec))
    
    % Carrier concentrations
    N = Con_vec(con_num);
    P = N;
    
    % Quasi Fermi levels
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
    E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);
    
    k_t_tag = R.k_t_mat; %[0:0.1:20]./L_z;
    theta_tag = [0:0.01:2*pi];
    
    z_exp = exp(-abs(repmat(R.z_grid,1,length(R.z_grid))-repmat(R.z_grid',length(R.z_grid),1)));
     
    for (ee = 1:length(E_exc))
        w = (E_exc(ee)./Consts.hbar);
        sum_TE = 0;
        sum_TM = 0;
        
        for (cb_num = 1:cb_len)
            for (vb_num = 1:vb_len)
                [cb_num, vb_num]
                
                % Carrier distributions -----------------------------------
                E_e = R.E_k_c{cb_num};
                E_h = R.E_k_v{vb_num};
                f_h = FermiDirac(E_fv, E_e, T);
                f_e = FermiDirac(E_fc, E_h, T);
                omega = f_e + f_h - 1;
                
                % Simulation definitions ----------------------------------
                mu_TE = sqrt(squeeze(R.rtsTE(vb_num,cb_num,:))).'.*Mb_s;
                mu_TM = sqrt(squeeze(R.rtsTM(vb_num,cb_num,:))).'.*Mb_s;
                Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
                Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
                Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
                Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
                Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
                
                for (kk = 1:length(R.k_t_mat))
                    for (kk_tag = 1:length(k_t_tag))
                        
                        % Calculate the q vector --------------------------
                        q_vec = sqrt(R.k_t_mat(kk)^2 + k_t_tag(kk_tag)^2 + 2*R.k_t_mat(kk)*k_t_tag(kk_tag).*cos(theta_tag));
                        q_vec(q_vec == 0) = 1e-5;
                        
                        % Theta matrix calculation ------------------------
                        C_Theta = Consts.e_0^2./(8*pi^2*R.Materials{2}.eps_r*Consts.eps_0);
                        
                        % Form factor for the various combinations
                        for (qq = 1:length(q_vec))
                            z_exp_q = z_exp.^q_vec(qq);
                            
                            % hole-electron
                            integrand_G_h = Xi_e.*Xi_h.*z_exp_q;
                            Gq_he(qq) = trapz(R.z_grid, trapz(R.z_grid,integrand_G_h,1), 2);
                            
                            % hole-hole
                            integrand_G_hh_h = Xi_h.*Xi_h.*z_exp_q;
                            Gq_hh(qq) = trapz(R.z_grid, trapz(R.z_grid,integrand_G_hh_h,1), 2);
                            
                            % electron-electron
                            integrand_G_ee = Xi_e.*Xi_e.*z_exp_q;
                            Gq_ee(qq) = trapz(R.z_grid, trapz(R.z_grid,integrand_G_ee,1), 2);
                        end
                        
                        % Lindhard formula for the dielectric function (for w=0)
                        for (qq = 1:length(q_vec))
                            for (kk = 1:length(R.k_t_mat))
                                theta = acos((R.k_t_mat(kk)^2 + q_vec.^2 - k_t_tag(kk_tag)^2)./(2*R.k_t_mat(kk).*q_vec));
                                for (tt = 1:length(theta))
                                    k_q_dist = sqrt(R.k_t_mat(kk).^2 + q_vec(qq)^2 - 2.*R.k_t_mat(kk)*q_vec(qq)*cos(theta(tt)));
                                    k_q_index = find(R.k_t_mat - k_q_dist < 1e-5);
                                    integrand_eps_theta_e(tt) = ((f_e(k_q_index)-f_e(kk))./(E_e(k_q_index)-E_e(kk)));
                                    integrand_eps_theta_h(tt) = ((f_h(k_q_index)-f_h(kk))./(E_h(k_q_index)-E_h(kk)));
                                end
                                integrand_eps_k_e(kk) = R.k_t_mat(kk)*trapz(theta, integrand_eps_theta_e);
                                integrand_eps_k_h(kk) = R.k_t_mat(kk)*trapz(theta, integrand_eps_theta_h);
                            end
                            eps_q(qq) = q_vec(qq)-((Consts.e_0^2)./(4*pi^2*R.Materials{2}.eps_r*Consts.eps_0)).*...
                                (trapz(R.k_t_mat, integrand_eps_k_e) + trapz(R.k_t_mat, integrand_eps_k_h));
                        end
                        
                        Theta(kk,kk_tag) = C_Theta.*trapz(theta_tag, Gq_eh_h./eps_q);
                        
                        % Coulomb-hole self energy -------------------------
                        C_E_CH = Consts.e_0^2./(4*pi^2*Materials{2}.eps_r*Consts.eps_0);
                        D_E_CH = C_E_CH.*trapz(q_vec, Gq_hh.*((1/(eps_q./q_vec))-1));
                        
                        % Screened-exchange shift energy -------------------
                        integrand_E_sx(kk_tag) = k_tag(kk_tag)*trapz(theta_tag,...
                            (1./eps_q).*(Gq_ee.*f_e(kk_tag)+Gq_hh.*f_h(kk_tag)));
                    end
                    
                    E_sx(kk) = -C_Theta*trapz(k_t_tag, integrand_E_sx);
                    
                    % Bandgap renormalizarion
                    w_tag(kk) = (E_h(kk) + E_e(kk) + D_E_CH + E_sx(kk))./Consts.hbar;
                end
                
                % Calculating the Coulomb enhancement factor --------------
                Xi_0_TE = (-1i/Consts.hbar).*((omega.*mu_TE)./(1i*(w_tag-w)+gamma));
                Xi_0_TM = (-1i/Consts.hbar).*((omega.*mu_TM)./(1i*(w_tag-w)+gamma));
                
                M_TE = eye(length(R.k_t_mat));
                M_TE = M_TE - repmat((dk./mu_TE).*Xi_0_TE, length(R.k_t_mat), 1).*Theta;
                M_TM = eye(length(R.k_t_mat));
                M_TM = M_TM - repmat((dk./mu_TM).*Xi_0_TM, length(R.k_t_mat), 1).*Theta;
                Q_k_TE{cb_num,vb_num} = M_TE\ones(length(R.k_t_mat));
                Q_k_TM{cb_num,vb_num} = M_TM\ones(length(R.k_t_mat));
                
                % Calculating gain sum elements ---------------------------
                I_TE = 2*trapz(R.k_t_mat, R.k_t_mat.*mu_TE^2.*omega.*(1./((w_mn-w).^2+gamma^2)).*Q_k_TE);
                I_TM = 2*trapz(R.k_t_mat, R.k_t_mat.*mu_TM^2.*omega.*(1./((w_mn-w).^2+gamma^2)).*Q_k_TM);
                sum_TE = sum_TE + I_TE;
                sum_TM = sum_TM + I_TM;
            end
        end
        
        G_TE(ii) = imag(C.*1i*w.*sum_TE*L_z);
        G_TM(ii) = imag(C.*1i*w.*sum_TM*L_z);
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

%                 Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
%                 Psi_h = abs(R.wf_h_mat{vb_num}(:,2));
%                 Psi_l = abs(R.wf_h_mat{vb_num}(:,3));
%                 for (ii=1:length(R.z_grid))
%                     int_h(ii) = trapz(R.z_grid, (Psi_e(ii).^2).*(Psi_h.^2).*exp(-q_vec(qq).*abs(R.z_grid(ii)-R.z_grid)));
%                     int_l(ii) = trapz(R.z_grid, (Psi_e(ii).^2).*(Psi_l.^2).*exp(-q_vec(qq).*abs(R.z_grid(ii)-R.z_grid)));
%                 end
%                 Gq_h(cb_num,vb_num,qq) = trapz(R.z_grid, int_h);
%                 Gq_l(cb_num,vb_num,qq) = trapz(R.z_grid, int_l);


%             % Screened potential calculation
%             Vq_h(cb_num,vb_num,:) = squeeze(Gq_h(cb_num,vb_num,:)).'.*(Consts.e_0^2./(2.*q_vec.*R.Materials{2}.eps_r));
%             Vq_l(cb_num,vb_num,:) = squeeze(Gq_l(cb_num,vb_num,:)).'.*(Consts.e_0^2./(2.*q_vec.*R.Materials{2}.eps_r));
%
%             figure(1);
%             subplot(211);
%             plot(q_vec, squeeze(Gq_h(cb_num,vb_num,:)), 'b', q_vec, squeeze(Gq_l(cb_num,vb_num,:)), 'r');
%             xlabel('q [m^{-1}]'); ylabel('G_q');
%             subplot(212);
%             plot(q_vec, squeeze(Vq_h(cb_num,vb_num,:)), 'b', q_vec, squeeze(Vq_l(cb_num,vb_num,:)), 'r');
%             xlabel('q [m^{-1}]'); ylabel('V_q');
%             drawnow;
%
%             % Lindhard formula calculation (for w=0)
%
%             for (qq = 1:length(q_vec))
%                 for (kk = 1:length(R.k_t_mat))
%                     for (tt = 1:length(theta))
%                         k_q_dist = sqrt(R.k_t_mat(kk).^2 + q_vec(qq)^2 - 2.*R.k_t_mat(kk)*q_vec(qq)*cos(theta(tt)));
%                         k_q_index = find(R.k_t_mat - k_q_dist(qq,kk,tt) < 1e-5);
%
%                         integrand_eps_theta_e(tt) = R.k_t_mat(kk)*((f_e(k_q_index)-f_e(kk))./(E_e(k_q_index)-E_e(kk)));
%                         integrand_eps_theta_h(tt) = R.k_t_mat(kk)*((f_h(k_q_index)-f_h(kk))./(E_h(k_q_index)-E_h(kk)));
%                     end
%                     integrand_eps_k_e(kk) = trapz(theta, integrand_eps_theta_e);
%                     integrand_eps_k_h(kk) = trapz(theta, integrand_eps_theta_h);
%                 end
%                 eps_q(qq) = 1-((Consts.e_0^2)./(4*pi^2*R.Materials{2}.eps_r*q_vec(qq))).*...
%                               (trapz(R.k_t_mat, integrand_eps_k_e) + trapz(R.k_t_mat, integrand_eps_k_h));
%             end

