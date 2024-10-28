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
x = 0.2;
T = 300;

Structure = { 'GaAlAs' , 100, x ;
              'GaAs', L_z*1e10, 0;
              'GaAlAs' , 100, x };

Res2 = GeneralStructure_TransferMatrix(Structure,T);

R = Res2;

close all;

%% Part 2 - Gain Spectra

% Definitions
colors = ['m';'c';'r';'g';'b';'y';'w';'k';'m';'c';'r';'g';'b';'y';'w';'k'];
              
con_vec = [2e10,3e10,5e10,7e10,2e11,3e11,5e11,7e11,3e12,5e12,7e12];    % [cm^-2]
gamma = 1e13;                        % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
%E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 0.001*Consts.e_0 : (R.Materials{2}.E_g+0.2)*Consts.e_0;
deltaE_vc_0 = (R.E_k_v{1}(1)/Consts.e_0 + R.E_k_c{1}(1)/Consts.e_0)
E_exc = (deltaE_vc_0-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (deltaE_vc_0+0.2)*Consts.e_0;
cb_len = 1; vb_len = 3;
init_E = 0.01*Consts.e_0;

% Simulation grids
sim_index = 1:300;
k_vec = R.k_t_mat(sim_index);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.1:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) -...
    repmat(R.z_grid',length(R.z_grid),1));

for (con_num = 1:length(con_vec))
    
    % Carrier concentrations
    N = con_vec(con_num);
    P = N;
    
    % Quasi Fermi levels
    %for (ii=1:length(R.E_k_v)), E_v_temp{ii} = -R.E_k_v{ii}; end
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
    E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);
    %E_fc = R.E_k_c{1}(1) + (R.E_k_c{1}(1) - E_fc);
    
    for (cb_num = 1:cb_len)
        E_e = R.E_k_c{cb_num}(sim_index);
        f_e = FermiDirac(E_fc, E_e, T);
        Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
        Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));        
        Xi_ee = Xi_e.*Xi_e.';
        Gq_ee = G_q(Xi_ee,k_vec,R.z_grid);
        G_Mat{cb_num,cb_num} = Gq_ee;
        
        for (vb_num = 1:vb_len)
            E_h = R.E_k_v{vb_num}(sim_index);
            f_h = FermiDirac(E_fv, E_h, T);
            Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
            Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));          
            Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);             
            Xi_hh = Xi_h.'.*Xi_h;
            Xi_eh = Xi_e.*Xi_h;
            
            Gq_eh = G_q(Xi_eh,k_vec,R.z_grid);
            Gq_hh = G_q(Xi_hh,k_vec,R.z_grid);         
            G_Mat{cb_num, cb_len+vb_num} = Gq_eh;
            G_Mat{cb_len+vb_num, cb_len+vb_num} = Gq_hh;
            
            eps_q = Eps_q(k_vec,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,R.Materials{2}.eps_r);
            eps_q(1) = eps_q(2);           
        end
    end
    
    k_min = k_vec(1);
    dk = k_vec(2)-k_vec(1);

    for (ii=1:cb_len)
        for (jj=1:vb_len)
            mu_TE(ii,jj,:) = sqrt(squeeze(R.rtsTE(jj,ii,sim_index)));
            mu_TM(ii,jj,:) = sqrt(squeeze(R.rtsTM(jj,ii,sim_index)));
        end
    end
        
    E_c = cell2mat(R.E_k_c(1:cb_len).'); 
    E_c = E_c(:,sim_index);
    E_v = cell2mat(R.E_k_v(1:vb_len).'); 
    %E_v = cell2mat(E_v_temp(1:vb_len).');
    E_v = E_v(:,sim_index);
    
    [G_TE_FCT, G_TM_FCT, R_sp_FCT, delta_n_FCT] = Gain_FCT(E_c, E_v, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, mu_TE, mu_TM, T, gamma, 'Lorentzian');
    [G_TE_MBT, G_TM_MBT, R_sp_MBT, delta_n_MBT, Delta_E_ch, Delta_E_sx_0, Q_mat] = Gain_MBT_HF(E_c, E_v, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, k_min, dk, k_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, G_Mat, eps_q, mu_TE, mu_TM, T, gamma);
    %[G_TE_MBT, G_TM_MBT, R_sp_MBT, delta_n_MBT, Delta_E_ch] = Gain_MBT_HF_Fast(E_c, E_v, E_fc, E_fv, cb_len, vb_len, k_min, dk, k_len, k_min, dk, k_len, k_min, dk, k_len, E_exc/Consts.hbar, length(E_exc), R.Materials, G_Mat, eps_q, mu_TE, mu_TM, T, gamma);

    Delta_E_ch_mat{con_num} = Delta_E_ch;
    Delta_E_sx_mat{con_num} = Delta_E_sx_0;
    
    index_TE = find(abs(G_TE_FCT) == max(abs(G_TE_FCT(1:round(length(G_TE_FCT))))));
    index_TM = find(abs(G_TM_FCT) == max(abs(G_TM_FCT(1:round(length(G_TM_FCT))))));
    index_Rsp = find(abs(R_sp_FCT) == max(abs(R_sp_FCT(1:round(length(R_sp_FCT))))));
    
    figure(2);
    subplot(311);
    plot(E_exc/Consts.e_0, G_TE_MBT, colors(con_num), E_exc/Consts.e_0, G_TE_FCT, strcat(colors(con_num),':')); hold on;
    plot(deltaE_vc_0*ones(size(G_TE_MBT)), G_TE_MBT, 'r');
    ylabel('G_T_E [cm^{-1}]'); grid on;
    %text(E_exc(index_TE)/Consts.e_0,G_TE_FCT(index_TE), [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    title(['x=' num2str(x) ', L_z=' num2str(L_z/1e-10) 'A, T=' num2str(T) 'K, [cb,vb]=[' num2str(cb_len) ',' num2str(vb_len) ']']);
    subplot(312);
    plot(E_exc/Consts.e_0, G_TM_MBT, colors(con_num), E_exc/Consts.e_0, G_TM_FCT, strcat(colors(con_num),':')); hold on;
    plot(deltaE_vc_0*ones(size(G_TM_MBT)), G_TM_MBT, 'r');
    %text(E_exc(index_TM)/Consts.e_0,G_TM_FCT(index_TM), [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    ylabel('G_T_M [cm^{-1}]'); grid on;
    subplot(313);
    plot(E_exc/Consts.e_0, R_sp_MBT, colors(con_num), E_exc/Consts.e_0, R_sp_FCT, strcat(colors(con_num),':')); hold on;
    plot(deltaE_vc_0*ones(size(R_sp_MBT)), R_sp_MBT, 'r');
    %text(E_exc(index_Rsp)/Consts.e_0,R_sp_FCT(index_Rsp), [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    xlabel('E [eV]'); ylabel('R_s_p'); grid on;
    drawnow;
end

