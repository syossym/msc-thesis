%function HartreeFock_Test3()

warning off;
global Consts;
global R G_TE G_TM z_exp;
format SHORT ENG;

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
con_vec = [2e9,4e9,6e9,2e10,4e10,6e10,8e10,1e11,2e11,3e11,4e11,6e11,8e11];    % [cm^-2]
%con_vec = 1e12;
gamma = 1e12;                                 % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
deltaE_vc_0 = (R.E_k_v{1}(1)/Consts.e_0 + R.E_k_c{1}(1)/Consts.e_0);
%E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 0.001*Consts.e_0 : (R.Materials{2}.E_g+0.4)*Consts.e_0;
E_exc = 1.45*Consts.e_0 : 1e-3*Consts.e_0 : 1.6*Consts.e_0;
E_g = GetMaterialBandGap(R.Materials{2}, T);   % [eV]
init_E = 0.001*Consts.e_0;
cb_len = 1; vb_len = 2;

% Simulation grids
sim_index = 1:50; % length(R.k_t_mat);
k_vec = interp(R.k_t_mat(sim_index), 3);
k_min = k_vec(1);
dk = k_vec(2)-k_vec(1);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.1:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) - repmat(R.z_grid',length(R.z_grid),1));

for (cb_num = 1:cb_len)
    for (vb_num = 1:vb_len)
        E_e_mat{cb_num} = interp1(R.k_t_mat(sim_index), R.E_k_c{cb_num}(sim_index), k_vec, 'pchip');
        E_h_mat{vb_num} = interp1(R.k_t_mat(sim_index), R.E_k_v{vb_num}(sim_index), k_vec, 'pchip');
        mu_TE_mat{cb_num,vb_num} = interp1(R.k_t_mat(sim_index), squeeze(R.rtsTE(vb_num,cb_num,sim_index)), k_vec, 'pchip');
        mu_TM_mat{cb_num,vb_num} = interp1(R.k_t_mat(sim_index), squeeze(R.rtsTM(vb_num,cb_num,sim_index)), k_vec, 'pchip');        
        E_e(cb_num,:) = E_e_mat{cb_num};
        E_h(vb_num,:) = E_h_mat{vb_num}; 
    end
end

disp('-- Calculating form factors');
eps_q = k_vec;
for (cb_num = 1:cb_len)
    Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
    Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
    Xi_ee = Xi_e.*Xi_e.';

    for (vb_num = 1:vb_len)
        disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
       
        Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
        Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
        Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
        Xi_hh = Xi_h.'.*Xi_h;
        Xi_eh = Xi_e.*Xi_h;
        
        Gq_ee{cb_num,vb_num} = ones(1,k_len); % G_q(Xi_ee,k_vec,R.z_grid);
        G_Mat{cb_num,cb_num} = Gq_ee{cb_num,vb_num};
        Gq_eh{cb_num,vb_num} = ones(1,k_len); % G_q(Xi_eh,k_vec,R.z_grid);
        Gq_hh{cb_num,vb_num} = ones(1,k_len); % G_q(Xi_hh,k_vec,R.z_grid);
        G_Mat{cb_num, cb_len+vb_num} = Gq_eh{cb_num,vb_num};
        G_Mat{cb_len+vb_num, cb_len+vb_num} = Gq_hh{cb_num,vb_num};
              
    end
end

for (con_num = 1:length(con_vec))
    
    % Carrier concentrations
    P = 1e6;
    N = con_vec(con_num) + P;
    
    disp(['-- Concentration - N=' num2str(N, '%1.0e') ' cm^-2, P=' num2str(P, '%1.0e') ' cm^-2']);
    
    %E_fc = 3.6e-11*N*1e-3*Consts.e_0 + R.E_k_c{1}(1);
    
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000);
    E_fv = fzero(@(E_fv) TargetFermiIntegral(k_vec, E_h_mat, E_fv, P*1e4, T, L_z, 'h'), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(k_vec, E_e_mat, E_fc, N*1e4, T, L_z, 'c'), init_E, options);
    
    E_f_vec(con_num) = E_fc;
    
    PlotQuasiFermiLevel(h_Band_Diagrams, E_fc, E_fv, k_vec, T, E_e_mat, E_h_mat, N, P);
    %PlotFermiLevel(h_Band_Diagrams, E_fc,  R.k_t_mat, T, R.E_k_c, R.E_k_v, N);
    
    close all;
    disp(['-- Parameter calculations']);
    disp('Calculating screened dielectric function');
    [eps_q, eps_q_fit] = CalculateScreening(k_vec, theta, E_fv, E_fc, T, E_e_mat(1:cb_len),E_h_mat(1:vb_len), R.Materials{2}.eps_r, Gq_ee, Gq_hh);
    
    for (cb_num = 1:cb_len)
        f_e = FermiDirac(E_fc, E_e_mat{cb_num}, T);
        for (vb_num = 1:vb_len)
            disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
            f_h = FermiDirac(E_fv, E_h_mat{vb_num}, T);
            disp('Calculating Theta matrix and Screened-exchange shift energy');
            [Theta{cb_num,vb_num}, D_E_SX{cb_num,vb_num}] = CalculateThetaESX(Gq_ee{cb_num,vb_num},Gq_hh{cb_num,vb_num},Gq_eh{cb_num,vb_num},theta,eps_q_fit,k_vec,f_e,f_h,R.Materials{2}.eps_r);
        end
    end
    
    disp('-- Spectrum calculation');   
    [G_TE_FCT, G_TM_FCT, R_sp_FCT, delta_n_FCT] = Gain_FCT(E_e, E_h, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, mu_TE_mat, mu_TM_mat, T, gamma, 'Lorentzian');
    [G_TE_MBT, G_TM_MBT, R_sp_MBT, delta_n_MBT, Delta_E_ch, Delta_E_sx_0, Q_mat, delta_renorm] = Gain_MBT_HF(E_e, E_h, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, k_min, dk, k_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, G_Mat, eps_q, mu_TE_mat, mu_TM_mat, T, gamma, Theta);
    
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

PlotAbsorptionSpEmission(E_exc, alpha_TE_FCT_vec, alpha_TM_FCT_vec, alpha_TE_MBT_vec, alpha_TM_MBT_vec, Rsp_FCT_vec, Rsp_MBT_vec, con_vec, T, E_f_vec, E_e_mat, E_h_mat, d_E, E_g);
