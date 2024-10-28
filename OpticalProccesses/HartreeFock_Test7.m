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
T = 5;
Structure = { 'GaAlAs' , 100, x ;
              'GaAs', L_z*1e10, 0;
              'GaAlAs' , 100, x };

%R = GeneralStructure_TransferMatrix(Structure, T);
%close all;
h_Band_Diagrams = PlotBandDiagram(R.k_t_mat, R.E_k_c, R.E_k_v, x, T, L_z, R.Materials{2});


%% Part 2 - Gain Spectra

% Definitions
con_vec = [2e10,6e10,2e11,6e11,2e12,4e12];    % [cm^-2]
%con_vec = 2e12;
gamma = 8e12;                                      % [s^-1]
dk = R.k_t_mat(2) - R.k_t_mat(1);
deltaE_vc_0 = (R.E_k_v{1}(1)/Consts.e_0 + R.E_k_c{1}(1)/Consts.e_0);
%E_exc =    (deltaE_vc_0-0.01)*Consts.e_0 : 5e-5*Consts.e_0 : (deltaE_vc_0+0.01)*Consts.e_0;
%E_exc = 1.51*Consts.e_0 : 5e-5*Consts.e_0 : 1.54*Consts.e_0;
E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 0.001*Consts.e_0 : (R.Materials{2}.E_g+0.3)*Consts.e_0;
Mb = sqrt(Consts.e_0*R.Materials{2}.E_p*Consts.m_0/6);
C_g = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
init_E = 0.01*Consts.e_0;
cb_len = 2; vb_len = 3;
h_Gain_Sp = figure('Name','Gain and Spontanious Emission');

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
    for (vb_num = 1:vb_len)  
        disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
        
        Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
        Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
        Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
        Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
        Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
        Xi_ee = Xi_e.*Xi_e.';
        Xi_hh = Xi_h.'.*Xi_h;
        Xi_eh = Xi_e.*Xi_h;
        
        Gq_eh{cb_num,vb_num} = G_q(Xi_eh,k_vec,R.z_grid);
        Gq_ee{cb_num,vb_num} = G_q(Xi_ee,k_vec,R.z_grid);
        Gq_hh{cb_num,vb_num} = G_q(Xi_hh,k_vec,R.z_grid);
    end
end

for (con_num = 1:length(con_vec))
    
    % Carrier concentrations
    N = con_vec(con_num);
    P = N;
    
    disp(['-- Concentration - N=' num2str(N, '%1.0e') ' cm^-2, P=' num2str(P, '%1.0e') ' cm^-2']);
    
    % Quasi Fermi levels
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000);
    E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z, 'c'), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z, 'h'), init_E, options);
    
    PlotQuasiFermiLevel(h_Band_Diagrams, E_fc, E_fv, R.k_t_mat, T, R.E_k_c, R.E_k_v, N, P);
    
    D_E_CH_con = {};
    D_E_SX_con = {};
    w_tag_con  = {};
    
    disp(['-- Parameter calculations']);
    for (cb_num = 1:cb_len)
        for (vb_num = 1:vb_len)
            disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
            
            E_e = R.E_k_c{cb_num}(sim_index);
            E_h = R.E_k_v{vb_num}(sim_index);
            f_h = FermiDirac(E_fv, E_h, T);
            f_e = FermiDirac(E_fc, E_e, T);

            disp('Calculating screened dielectric function');
            eps_q = Eps_q(k_vec,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,R.Materials{2}.eps_r);
            eps_q(1) = eps_q(2);
            
            disp('Calculating Theta matrix and Screened-exchange shift energy');
            [Theta{cb_num,vb_num}, D_E_SX{cb_num,vb_num}] = ThetaESxCalculation(Gq_ee{cb_num,vb_num},Gq_hh{cb_num,vb_num},Gq_eh{cb_num,vb_num},theta,eps_q,k_vec,f_e,f_h,R.Materials{2}.eps_r);
            
            disp('Calculating D_E_CH');
            q_vec = k_vec; q_vec(1) = q_vec(2);
            D_E_CH{cb_num,vb_num} = E_CH(Gq_hh{cb_num,vb_num},eps_q,k_vec,q_vec,R.Materials{2}.eps_r);
                     
            disp('Calculating bandgaps');
            w_mn_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num})./Consts.hbar;
            w_mn{cb_num,vb_num} = (E_h + E_e)./Consts.hbar;
        end
    end
    
    disp('-- Spectrum calculation');
    for (ee = 1:length(E_exc))
        E_exc(ee);
        w = (E_exc(ee)./Consts.hbar);
        sum_TE_MBT = 0;
        sum_TM_MBT = 0;
        sum_TE_FCT = 0;
        sum_TM_FCT = 0;
        sum_Rsp_FCT = 0;
        sum_Rsp_MBT = 0;
        
        for (cb_num = 1:cb_len)
            for (vb_num = 1:vb_len)
              
                %disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
                
                % Carrier distributions -----------------------------------
                E_e = R.E_k_c{cb_num}(sim_index);
                E_h = R.E_k_v{vb_num}(sim_index);
                f_h = FermiDirac(E_fv, E_h, T);
                f_e = FermiDirac(E_fc, E_e, T);
                omega = f_e + f_h - 1;
                
                % Simulation definitions ----------------------------------
                mu_TE = sqrt(squeeze(R.rtsTE(vb_num,cb_num,sim_index))).'.*Mb;
                mu_TM = sqrt(squeeze(R.rtsTM(vb_num,cb_num,sim_index))).'.*Mb;
                
                %disp('Calculating Coulomb enhancement factor');
                Xi_0_TE = (-1i/Consts.hbar).*((omega.*mu_TE)./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma));
                Xi_0_TM = (-1i/Consts.hbar).*((omega.*mu_TM)./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma));
                M_TE = repmat((1./mu_TE).', 1, k_len).*repmat(Xi_0_TE.*k_vec, k_len, 1).*Theta{cb_num,vb_num};
                Q_k_TE = 1./(1-trapz(k_vec, M_TE, 2));
                M_TM = repmat((1./mu_TM).', 1, k_len).*repmat(Xi_0_TM.*k_vec, k_len, 1).*Theta{cb_num,vb_num};
                Q_k_TM = 1./(1-trapz(k_vec, M_TM, 2));
                Q_k_Rsp = 1./(1-(2/3)*trapz(k_vec, M_TE, 2)+(1/3)*trapz(k_vec, M_TM, 2));

                %disp('Calculating gain sum elements');
                I_TE_FCT = 2*trapz(k_vec, k_vec.*mu_TE.^2.*(omega).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2)));
                I_TM_FCT = 2*trapz(k_vec, k_vec.*mu_TM.^2.*(omega).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2)));
                sum_TE_FCT = sum_TE_FCT + I_TE_FCT;
                sum_TM_FCT = sum_TM_FCT + I_TM_FCT;
                
                I_Rsp_FCT = 2*trapz(k_vec, k_vec.*((2/3)*mu_TE.^2+(1/3)*mu_TM.^2).*(f_e.*f_h).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2)));
                sum_Rsp_FCT = sum_Rsp_FCT + I_Rsp_FCT;
                
                I_TE_MBT = 2*trapz(k_vec, k_vec.*mu_TE.^2.*omega.*(1./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TE.');
                I_TM_MBT = 2*trapz(k_vec, k_vec.*mu_TM.^2.*omega.*(1./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TM.');
                sum_TE_MBT = sum_TE_MBT + I_TE_MBT;
                sum_TM_MBT = sum_TM_MBT + I_TM_MBT;
                
                I_Rsp_MBT = 2*trapz(k_vec, k_vec.*((2/3)*mu_TE.^2+(1/3)*mu_TM.^2).*(f_e.*f_h).*(1./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_Rsp.');
                sum_Rsp_MBT = sum_Rsp_MBT + I_Rsp_MBT;
            end
        end
        
        G_TE_FCT(ee) = C_g.*w.*sum_TE_FCT*L_z;
        G_TM_FCT(ee) = C_g.*w.*sum_TM_FCT*L_z;
        G_Rsp_FCT(ee) = C_g.*w.*sum_Rsp_FCT*L_z; 
        G_TE_MBT(ee) = C_g*L_z*imag(1i*w.*sum_TE_MBT);
        G_TM_MBT(ee) = C_g*L_z*imag(1i*w.*sum_TM_MBT);
        G_Rsp_MBT(ee) = C_g*L_z*imag(1i*w.*sum_Rsp_MBT);
        D_E_CH_con = {D_E_CH_con, D_E_CH};
        D_E_SX_con = {D_E_SX_con, D_E_SX};
        %w_tag_con = {w_tag_con, w_tag};
    end
    
    E_f = 3.6e-11*N/1000;     % eV
    E_g = GetMaterialBandGap(R.Materials{2}, T); % eV
    E_f = E_f + (R.E_k_c{1}(1)/Consts.e_0);
    
    index_TE_MBT = find(abs(G_TE_MBT) == max(abs(G_TE_MBT(1:round(length(G_TE_MBT)/3)))));
    index_TM_MBT = find(abs(G_TM_MBT) == max(abs(G_TM_MBT(1:round(length(G_TM_MBT)/3)))));
    index_Rsp_MBT = find(abs(G_Rsp_MBT) == max(abs(G_Rsp_MBT(1:round(length(G_Rsp_MBT)/3)))));
       
    figure(h_Gain_Sp);
    subplot(311);
    plot(E_exc/Consts.e_0, G_TE_MBT/100, 'b', E_exc/Consts.e_0, G_TE_FCT/100, 'r:'); hold on;
    %plot(deltaE_vc_0*ones(size(G_TE_MBT)), G_TE_MBT/100, 'g:');
    plot(deltaE_vc_0*ones(size(G_TE_MBT)), G_TE_MBT/100, 'g:', E_f*ones(size(G_TE_MBT)), G_TE_MBT/100, 'y', ...
         E_g*ones(size(G_TE_MBT)), G_TE_MBT/100, 'c', (R.E_k_c{1}(1)/Consts.e_0)*ones(size(G_TE_MBT)), G_TE_MBT/100, 'k', ...
         (R.E_k_c{2}(1)/Consts.e_0)*ones(size(G_TE_MBT)), G_TE_MBT/100, 'm');   
    ylabel('G_{TE} [cm^{-1}]');
    text(E_exc(index_TE_MBT)/Consts.e_0,G_TE_MBT(index_TE_MBT)/100, [num2str(N,'%1.0e'), ' cm^{-2}'], 'FontSize', 8);
    subplot(312);
    plot(E_exc/Consts.e_0, G_TM_MBT/100, 'b', E_exc/Consts.e_0, G_TM_FCT/100, 'r:'); hold on;
    %plot(deltaE_vc_0*ones(size(G_TM_MBT)), G_TM_MBT/100, 'g:');
    plot(deltaE_vc_0*ones(size(G_Rsp_MBT)), G_Rsp_MBT/100, 'g:', E_f*ones(size(G_Rsp_MBT)), G_TM_MBT/100, 'y', ...
         E_g*ones(size(G_TM_MBT)), G_TM_MBT/100, 'c', (R.E_k_c{1}(1)/Consts.e_0)*ones(size(G_TM_MBT)), G_TM_MBT/100, 'k', ...
         (R.E_k_c{2}(1)/Consts.e_0)*ones(size(G_TM_MBT)), G_TM_MBT/100, 'm');   
    text(E_exc(index_TM_MBT)/Consts.e_0,G_TM_MBT(index_TM_MBT)/100, [num2str(N,'%1.0e'), ' cm^{-2}'], 'FontSize', 8);
    ylabel('G_{TM} [cm^{-1}]');
    subplot(313);
    plot(E_exc/Consts.e_0, G_Rsp_MBT, 'b', E_exc/Consts.e_0, G_Rsp_FCT, 'r:'); hold on;
    plot(deltaE_vc_0*ones(size(G_Rsp_MBT)), G_Rsp_MBT, 'g:', E_f*ones(size(G_Rsp_MBT)), G_Rsp_MBT, 'y', ...
         E_g*ones(size(G_Rsp_MBT)), G_Rsp_MBT, 'c', (R.E_k_c{1}(1)/Consts.e_0)*ones(size(G_Rsp_MBT)), G_Rsp_MBT, 'k', ...
         (R.E_k_c{2}(1)/Consts.e_0)*ones(size(G_Rsp_MBT)), G_Rsp_MBT, 'm');   
    text(E_exc(index_Rsp_MBT)/Consts.e_0, G_Rsp_MBT(index_Rsp_MBT), [num2str(N,'%1.0e'), ' cm^{-2}'], 'FontSize', 8);
    ylabel('R_{sp} [a.u]');
    xlabel('E [eV]'); ylabel('R_{sp} [a.u.]');
    drawnow;
end
