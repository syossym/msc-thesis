%function HartreeFock_Test4()

warning off;
global Consts;
global R G_TE G_TM z_exp;

%% Part 1 - Subband energies
% close all; clc; %clear all;
% warning off;
%
% %Init project
global project_path;
project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
cd(project_path);
run('.\Common\AddPath.m');

Constants;

% well_width = 200   % [A]
% x = 1;
% Res1 = SingleQW_kp(well_width,x);

L_z = 200e-10;

Structure = { 'GaAlAs' , 100, 0.2 ;
    'GaAs', L_z*1e10, 0;
    'GaAlAs' , 100, 0.2 };

Res2 = GeneralStructure_TransferMatrix(Structure);

R = Res2;

close all;

%% Part 2 - Gain Spectra

% Definitions
con_vec = [2e11,6e11,2e12,3e12,4e12,6e12,8e12];    % [cm^-2]
T = 10;
gamma = 1e13;                                      % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 0.001*Consts.e_0 : (R.Materials{2}.E_g+0.5)*Consts.e_0;
Mb_s = Consts.e_0*R.Materials{2}.E_p*Consts.m_0/6;
C_g = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
init_E = 0.01*Consts.e_0;
cb_len = 1; vb_len = 3;

% Simulation grids
sim_index = 1:500;
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
    
    D_E_CH_con = {};
    D_E_SX_con = {};
    w_tag_con  = {};
    for (cb_num = 1:cb_len)
        for (vb_num = 1:vb_len)
            disp(['-- Parameter calculations - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
            
            % Carrier distributions -----------------------------------
            E_e = R.E_k_c{cb_num}(sim_index);
            E_h = R.E_k_v{vb_num}(sim_index);
            f_h = FermiDirac(E_fv, E_h, T);
            f_e = FermiDirac(E_fc, E_e, T);
            
            % Simulation definitions ----------------------------------
            Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
            Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
            Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
            Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
            Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
            Xi_ee = Xi_e.*Xi_e.';
            Xi_hh = Xi_h.'.*Xi_h;
            Xi_eh = Xi_e.*Xi_h;
            
            % Prelimenary clalculations -------------------------------
            disp('Calculating form factors');
            Gq_eh = G_q(Xi_eh,k_vec,R.z_grid);
            Gq_ee = G_q(Xi_ee,k_vec,R.z_grid);
            Gq_hh = G_q(Xi_hh,k_vec,R.z_grid);
            disp('Calculating screened dielectric function');
            eps_q = Eps_q(k_vec,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,R.Materials{2}.eps_r);
            eps_q(1) = eps_q(2);
            
            figure(1);
            subplot(211);
            plot(k_vec, Gq_eh, 'b', k_vec, Gq_ee, 'r', k_vec, Gq_hh, 'g');
            ylabel('G_q');
            legend('G_q^{eh}', 'G_q^{ee}', 'G_q^{hh}');
            subplot(212);
            plot(k_vec, eps_q./k_vec);
            xlabel('k'); ylabel('\epsilon_q');
            drawnow;
            
            % Theta matrix and Screened-exchange shift energy calculation --------------------------------
            disp('Calculating Theta matrix and Screened-exchange shift energy');
            [Theta{cb_num,vb_num}, D_E_SX{cb_num,vb_num}] = ThetaESxCalculation(Gq_ee,Gq_hh,Gq_eh,theta,eps_q,k_vec,f_e,f_h,R.Materials{2}.eps_r);
            
            % Xi matrix calculation -----------------------------------
            disp('Calculating D_E_CH');
            q_vec = k_vec; q_vec(1) = q_vec(2);
            D_E_CH{cb_num,vb_num} = E_CH(Gq_hh,eps_q,k_vec,q_vec,R.Materials{2}.eps_r);
            
            % Bandgap renormalization
            w_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num})./Consts.hbar;
        end
    end
    
    disp('-- Spectrum calculation');
    sum_TE = 0;
    sum_TM = 0;
    sum_sp = 0;
    sum_d_n = 0;
    for (cb_num = 1:cb_len)
        for (vb_num = 1:vb_len)
            
            % Carrier distributions -----------------------------------
            E_e = R.E_k_c{cb_num}(sim_index);
            E_h = R.E_k_v{vb_num}(sim_index);
            f_h = FermiDirac(E_fv, E_h, T);
            f_e = FermiDirac(E_fc, E_e, T);
            omega = f_e + f_h - 1;
            
            % Matrix elements -----------------------------------------
            mu_TE = sqrt(squeeze(R.rtsTE(vb_num,cb_num,sim_index)).*Mb_s).';
            mu_TM = sqrt(squeeze(R.rtsTM(vb_num,cb_num,sim_index)).*Mb_s).';
            
            for (ee = 1:length(E_exc))
                
                w = (E_exc(ee)./Consts.hbar);             
                C0 = C_g/w;
                
                [Qk_TE, Qk_TM, Qk_sp] = CalculateEnhancement(k_vec, theta, omega, D_E_CH, D_E_SX, w_tag, eps_q, mu_TE, mu_TM, gamma);
                
                g_TE(ee) = 2*C0*trapz(k_vec, real(Qk_TE));
                g_TM(ee) = 2*C0*trapz(k_vec, real(Qk_TM));
                r_sp(ee) = 2*C0*trapz(k_vec, real(Qk_sp));
                delta_n(ee) = (Consts.c*C0/w)*trapz(k_vec, -imag(Qk_TE));
                
                sum_TE = sum_TE + g_TE;
                sum_TM = sum_TM + g_TM;
                sum_sp = sum_sp + r_sp;
                sum_d_n = sum_d_n + delta_n;
                
                %                 Xi_0_TE = (-1i/Consts.hbar).*((omega.*mu_TE)./(1i*(w_tag{cb_num,vb_num}-w)+gamma));
                %                 Xi_0_TM = (-1i/Consts.hbar).*((omega.*mu_TM)./(1i*(w_tag{cb_num,vb_num}-w)+gamma));
                %
                %                 % Coulomb enhancement factor calculation ------------------
                %                 M_TE = eye(k_len);
                %                 M_TE = M_TE - repmat((dk./mu_TE).', 1, k_len).*repmat(Xi_0_TE, k_len, 1).*Theta{cb_num,vb_num}*1e9;%(:,end:-1:1).*1e8;
                %                 M_TM = eye(k_len);
                %                 M_TM = M_TM - repmat((dk./mu_TM).', 1, k_len).*repmat(Xi_0_TM, k_len, 1).*Theta{cb_num,vb_num}*1e9;%(:,end:-1:1).*1e8;
                %                 Q_k_TE = M_TE\ones(k_len,1);
                %                 Q_k_TM = M_TM\ones(k_len,1);
                %
                %                 % Calculating gain sum elements ---------------------------
                %                 I_TE = 2*trapz(k_vec, k_vec.*mu_TE.^2.*omega.*(1./(1i*(w_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TE.');
                %                 I_TM = 2*trapz(k_vec, k_vec.*mu_TM.^2.*omega.*(1./(1i*(w_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TM.');
                
            end
        end
        
        %G_TE(ee) = imag(C_g.*1i*w.*sum_TE*L_z);
        %G_TM(ee) = imag(C_g.*1i*w.*sum_TM*L_z);
        
    end
    
    index_TE = find(abs(G_TE) == max(abs(G_TE(1:round(length(G_TE)/3)))));
    index_TM = find(abs(G_TM) == max(abs(G_TM(1:round(length(G_TM)/3)))));
    
    figure(2);
    subplot(211);
    plot(E_exc/Consts.e_0, G_TE/100); hold on;
    plot(R.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
    ylabel('G [cm^{-1}]');
    text(E_exc(index_TE)/Consts.e_0,G_TE(index_TE)/100, [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    subplot(212);
    plot(E_exc/Consts.e_0, G_TM/100); hold on;
    plot(R.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
    text(E_exc(index_TM)/Consts.e_0,G_TM(index_TM)/100, [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    xlabel('E [eV]'); ylabel('G [cm^{-1}]');
    
    D_E_CH_con = {D_E_CH_con, D_E_CH};
    D_E_SX_con = {D_E_SX_con, D_E_SX};
    w_tag_con = {w_tag_con, w_tag};
end
