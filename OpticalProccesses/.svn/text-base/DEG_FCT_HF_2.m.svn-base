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
x = 0.1;
T = 2;

Structure = { 'GaAlAs' , 100, x ;
              'GaAs', L_z*1e10, 0;
              'GaAlAs' , 100, x };
          
%R = GeneralStructure_TransferMatrix(Structure, T);
close all;
h_Band_Diagrams = PlotBandDiagram(R.k_t_mat, R.E_k_c, R.E_k_v, x, T, L_z, R.Materials{2});

%% Part 2 - Gain Spectra

% Definitions
colors = ['m';'c';'r';'g';'b';'y';'w';'k';'m';'c';'r';'g';'b';'y';'w';'k'];

con_vec = [2e10,4e10,6e10,2e11,4e11,6e11,4e12,6e12];    % [cm^-2]
gamma = 1e12;                        % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
%E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 0.001*Consts.e_0 : (R.Materials{2}.E_g+0.2)*Consts.e_0;
deltaE_vc_0 = (R.E_k_v{1}(1)/Consts.e_0 + R.E_k_c{1}(1)/Consts.e_0);
%E_exc = (deltaE_vc_0-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (deltaE_vc_0+0.2)*Consts.e_0;
E_exc = 1.4*Consts.e_0 : 1e-4*Consts.e_0 : 1.6*Consts.e_0;
E_g = GetMaterialBandGap(R.Materials{2}, T); % eV
cb_len = 1; vb_len = 1;
init_E = 0.01*Consts.e_0;

% Simulation grids
sim_index = 1:length(R.k_t_mat);
k_vec = R.k_t_mat(sim_index);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.1:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) -...
    repmat(R.z_grid',length(R.z_grid),1));

disp('-- Calculating form factors');
for (cb_num = 1:cb_len)
    E_e = R.E_k_c{cb_num}(sim_index);
    Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
    Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
    Xi_ee = Xi_e.*Xi_e.';
    Gq_ee = ones(1, k_len); % G_q(Xi_ee,k_vec,R.z_grid);
    G_Mat{cb_num,cb_num} = Gq_ee;
    
    for (vb_num = 1:vb_len)
        disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
        
        E_h = R.E_k_v{vb_num}(sim_index);
        Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
        Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
        Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
        Xi_hh = Xi_h.'.*Xi_h;
        Xi_eh = Xi_e.*Xi_h;
        
        Gq_eh = G_q(Xi_eh,k_vec,R.z_grid);
        Gq_hh = G_q(Xi_hh,k_vec,R.z_grid);
        G_Mat{cb_num, cb_len+vb_num} = ones(1, k_len); % Gq_eh;
        G_Mat{cb_len+vb_num, cb_len+vb_num} = ones(1, k_len); % Gq_hh;
        
    end
end

for (con_num = 1:length(con_vec))
    
    % Carrier concentrations
    P = 1e6;
    N = con_vec(con_num) + P;
    
    disp(['-- Concentration - N=' num2str(N, '%1.0e') ' cm^-2, P=' num2str(P, '%1.0e') ' cm^-2']);
    
    % Quasi Fermi levels
    %E_fc = 3.6e-11*N*1e-3*Consts.e_0 + R.E_k_c{1}(1);
    
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000);
    E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z, 'h'), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z, 'c'), init_E, options);
    E_f_vec(con_num) = E_fc;
       
    PlotQuasiFermiLevel(h_Band_Diagrams, E_fc, E_fv, R.k_t_mat, T, R.E_k_c, R.E_k_v, N, P);
    
    k_vec_interp = k_vec; % [min(k_vec) : (k_vec(2)-k_vec(1)) : max(k_vec)];
    for (cb_num = 1:cb_len)
        f_e = FermiDirac(E_fc, E_e, T);
        E_e = R.E_k_c{cb_num}(sim_index);
        E_e_mat{cb_num} = E_e;
        
        for (vb_num = 1:vb_len)           
            disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
            
            E_h = R.E_k_v{vb_num}(sim_index);
            f_h = FermiDirac(E_fv, E_h, T);
            E_h_mat{vb_num} = E_h;
            
            disp(['  Calculating screened dielectric function']);
            E_e = interp1(k_vec, E_e, k_vec_interp, 'pchip');
            E_h = interp1(k_vec, E_h, k_vec_interp, 'pchip');
            f_e = interp1(k_vec, f_e, k_vec_interp, 'pchip');
            f_h = interp1(k_vec, f_h, k_vec_interp, 'pchip');
            [eps_q, eps_q_fit] = CalculateScreening(k_vec, theta, E_fv, E_fc, T, E_e_mat(1:cb_len),E_h_mat(1:vb_len), R.Materials{2}.eps_r, Gq_ee, Gq_hh);

            %eps_q = Eps_q(k_vec_interp,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,R.Materials{2}.eps_r);
            %eps_q = k_vec_interp;
            eps_q(1) = eps_q(2);
        end
    end
   
    k_min = k_vec_interp(1);
    dk = k_vec_interp(2)-k_vec_interp(1);
    k_len = length(k_vec_interp);
    
    for (ii=1:cb_len)
        for (jj=1:vb_len)
            mu_TE_temp = sqrt(squeeze(R.rtsTE(jj,ii,sim_index))).';
            mu_TM_temp = sqrt(squeeze(R.rtsTM(jj,ii,sim_index))).';
            mu_TE(ii,jj,:) = interp1(k_vec, mu_TE_temp, k_vec_interp, 'pchip');
            mu_TM(ii,jj,:) = interp1(k_vec, mu_TM_temp, k_vec_interp, 'pchip');
        end
    end
    
    E_c = cell2mat(R.E_k_c(1:cb_len).');
    E_c = E_c(:,sim_index);
    E_c = interp1(k_vec, E_c, k_vec_interp);
    E_v = cell2mat(R.E_k_v(1:vb_len).');
    %E_v = cell2mat(E_v_temp(1:vb_len).');
    E_v = E_v(:,sim_index);
    E_v = interp1(k_vec, E_v, k_vec_interp);
    
    disp('-- Spectrum calculation');
    [G_TE_FCT, G_TM_FCT, R_sp_FCT, delta_n_FCT] = Gain_FCT(E_c, E_v, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, mu_TE, mu_TM, T, gamma, 'Lorentzian');
    [G_TE_MBT, G_TM_MBT, R_sp_MBT, delta_n_MBT, Delta_E_ch, Delta_E_sx_0, Q_mat, delta_renorm] = Gain_MBT_HF(E_c, E_v, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, k_min, dk, k_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, G_Mat, eps_q, mu_TE, mu_TM, T, gamma);
    
    alpha_TE_FCT_vec{con_num} = -G_TE_FCT;
    alpha_TM_FCT_vec{con_num} = -G_TM_FCT;
    alpha_TE_MBT_vec{con_num} = -G_TE_MBT;
    alpha_TM_MBT_vec{con_num} = -G_TM_MBT;
    Rsp_FCT_vec{con_num} = R_sp_FCT;
    Rsp_MBT_vec{con_num} = R_sp_MBT;
    delta_n_FCT_vec{con_num} = delta_n_FCT;
    delta_n_MBT_vec{con_num} = delta_n_MBT;
    
    Delta_E_ch_mat{con_num} = Delta_E_ch;
    Delta_E_sx_mat{con_num} = Delta_E_sx_0;
    d_E(con_num) = delta_renorm{1,1}(1);
end

PlotAbsorptionSpEmission(E_exc, alpha_TE_FCT_vec, alpha_TM_FCT_vec, alpha_TE_MBT_vec, alpha_TM_MBT_vec, Rsp_FCT_vec, Rsp_MBT_vec, con_vec, T, E_f_vec, R.E_k_c, R.E_k_v, d_E, E_g);