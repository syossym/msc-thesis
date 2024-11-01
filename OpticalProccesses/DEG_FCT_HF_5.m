
warning off; clc;

%% Init Simulation

global Consts;
global R G_TE G_TM z_exp;
global project_path;

project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
cd(project_path);
run('.\Common\AddPath.m');

Constants;

%% Structure Simulation

% Structure defintion
L_z = 200e-10;
x = 0.1;
T = 2;
Structure = { 'GaAlAs' , 100, x ;
              'GaAs', L_z*1e10, 0;
              'GaAlAs' , 100, x };

%R = GeneralStructure_TransferMatrix(Structure, T);
close all;
h_Band_Diagrams = PlotBandDiagram(R.k_t_mat, R.E_k_c, R.E_k_v, x, T, L_z, R.Materials{2});

%% Optical Simulation

% Definitions -------------------------------------------------------------
con_vec = [2e9,4e9,6e9,8e9,2e10,4e10,6e10,8e10,2e11,4e11,6e11,8e11];    % [cm^-2]
gamma   = 7e11;                                                         % [s^-1]
Mb      = sqrt(Consts.e_0*R.Materials{2}.E_p*Consts.m_0/6);
E_p     = Consts.e_0*R.Materials{2}.E_p;                                % [J] 
E_g     = GetMaterialBandGap(R.Materials{2}, T);                        % [eV]
init_E  = 0.001*Consts.e_0;

C_g = 1/(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
C_d_n = 1/(Consts.hbar*Consts.eps_0*L_z*pi*R.Materials{2}.n);       
C_rsp = R.Materials{2}.n/(pi^3*Consts.hbar^2*Consts.c^3*L_z*Consts.eps_0);

deltaE_vc_0 = (R.E_k_v{1}(1)/Consts.e_0 + R.E_k_c{1}(1)/Consts.e_0);
E_exc       = 1.4*Consts.e_0 : 1e-4*Consts.e_0 : 1.75*Consts.e_0;
               
cb_len = 1; 
vb_len = 1;

% Simulation grids --------------------------------------------------------
sim_index = 1:50; %length(R.k_t_mat);
dk = R.k_t_mat(2) - R.k_t_mat(1);
k_vec = interp(R.k_t_mat(sim_index), 2);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.01:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) - repmat(R.z_grid',length(R.z_grid),1));

% Variable pre-processing -------------------------------------------------
for (cb_num = 1:cb_len)
    for (vb_num = 1:vb_len)
        E_e_mat{cb_num} = interp1(R.k_t_mat(sim_index), R.E_k_c{cb_num}(sim_index), k_vec, 'pchip');
        E_h_mat{vb_num} = interp1(R.k_t_mat(sim_index), R.E_k_v{vb_num}(sim_index), k_vec, 'pchip');
        mu_TE_mat{cb_num,vb_num} = interp1(R.k_t_mat(sim_index), squeeze(R.rtsTE(vb_num,cb_num,sim_index)), k_vec, 'pchip');
        mu_TM_mat{cb_num,vb_num} = interp1(R.k_t_mat(sim_index), squeeze(R.rtsTM(vb_num,cb_num,sim_index)), k_vec, 'pchip');
    end
end

% Form factor calculation -------------------------------------------------
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

% Main simulation ---------------------------------------------------------
for (con_num = 1:length(con_vec))
    t_start = tic;
    
    % Carrier concentrations
    P = 1e6;                      % [cm^-2]
    N = con_vec(con_num) + P;
    
    disp(['-- Concentration - N=' num2str(N, '%1.0e') ' cm^-2, P=' num2str(P, '%1.0e') ' cm^-2']);
    
    %E_fc = 3.6e-11*N*1e-3*Consts.e_0 + R.E_k_c{1}(1);
    
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000);
    E_fv = fzero(@(E_fv) TargetFermiIntegral(k_vec, E_h_mat, E_fv, P*1e4, T, L_z, 'h'), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(k_vec, E_e_mat, E_fc, N*1e4, T, L_z, 'c'), init_E, options);
    
    E_f_vec(con_num) = E_fc;
    
    PlotQuasiFermiLevel(h_Band_Diagrams, E_fc, E_fv, k_vec, T, E_e_mat, E_h_mat, N, P);
    %PlotFermiLevel(h_Band_Diagrams, E_fc,  R.k_t_mat, T, R.E_k_c, R.E_k_v, N);
    
    D_E_CH_con = {};
    D_E_SX_con = {};
    w_tag_con  = {};

    disp(['-- Parameter calculations']);
    disp('Calculating screened dielectric function');
    [eps_q, eps_q_fit] = CalculateScreening(k_vec, theta, E_fv, E_fc, T, E_e_mat(1:cb_len),E_h_mat(1:vb_len), R.Materials{2}.eps_r, Gq_ee, Gq_hh);
    
    for (cb_num = 1:cb_len)
        for (vb_num = 1:vb_len)
            disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);
            
            E_e = E_e_mat{cb_num};
            E_h = E_h_mat{vb_num};
            
            %f_h4 = 1 - exp(-(E_h-E_fc)/(Consts.k_B*T));
            f_h = FermiDirac(E_fv, E_h, T);
            f_e = FermiDirac(E_fc, E_e, T);
            
            disp('Calculating Theta matrix and Screened-exchange shift energy');
            [Theta{cb_num,vb_num}, D_E_SX{cb_num,vb_num}] = CalculateThetaESX(Gq_ee{cb_num,vb_num},Gq_hh{cb_num,vb_num},Gq_eh{cb_num,vb_num},theta,eps_q_fit,k_vec,f_e,f_h,R.Materials{2}.eps_r);
            
            disp('Calculating D_E_CH');
            D_E_CH{cb_num,vb_num} = CalculateECH(Gq_hh{cb_num,vb_num},eps_q_fit,k_vec,R.Materials{2}.eps_r);
            
            disp('Calculating bandgaps');
            w_mn_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num})./Consts.hbar;
            %w_mn_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)))./Consts.hbar;
            w_mn{cb_num,vb_num} = (E_h + E_e)./Consts.hbar;
            
            delta_renorm{cb_num,vb_num} = (D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num});
            %delta_renorm{cb_num,vb_num} = (D_E_CH{cb_num,vb_num}.*ones(size(k_vec)));
            d_E(con_num,cb_num,vb_num) = delta_renorm{cb_num,vb_num}(1);
        end
    end
    
    disp('-- Spectrum calculation');
    %h_FCT_Integrands = figure;
    %h_HF_Integrands = figure;
    for (ee = 1:length(E_exc))
        w = (E_exc(ee)./Consts.hbar);
        mu_0 = (Consts.e_0/(sqrt(3)*w))*sqrt(E_p/2/Consts.m_0);
        
        sum_TE_HF = 0;
        sum_TM_HF = 0;
        sum_TE_N_HF = 0;
        sum_TM_N_HF = 0;
        sum_TE_FCT = 0;
        sum_TM_FCT = 0;
        sum_TE_N_FCT = 0;
        sum_TM_N_FCT = 0;
        sum_Rsp_FCT = 0;
        sum_Rsp_TE_FCT = 0; 
        sum_Rsp_TM_FCT = 0;
        sum_Rsp_HF = 0;
        sum_Rsp_TE_HF = 0;
        sum_Rsp_TM_HF = 0;
        
        for (cb_num = 1:cb_len)
            for (vb_num = 1:vb_len)
                
                % Carrier distributions 
                E_e = E_e_mat{cb_num};
                E_h = E_h_mat{vb_num};
                f_h = FermiDirac(E_fv, E_h, T);
                f_e = FermiDirac(E_fc, E_e, T);
                omega_g = f_e + f_h - 1;
                omega_rsp = f_e.*f_h;
                
                % Simulation definitions
                mu_TE = sqrt(mu_TE_mat{cb_num,vb_num}).*mu_0;
                mu_TM = sqrt(mu_TM_mat{cb_num,vb_num}).*mu_0;
                   
                % Free carrier theory 
                int_FCT_1 = k_vec.*mu_TE.^2.*(omega_g).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_2 = k_vec.*mu_TM.^2.*(omega_g).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_3 = k_vec.*mu_TE.^2.*(omega_g).*((w_mn{cb_num,vb_num}-w)./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_4 = k_vec.*mu_TM.^2.*(omega_g).*((w_mn{cb_num,vb_num}-w)./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_5 = k_vec.*((2/3)*mu_TE.^2+(1/3)*mu_TM.^2).*(omega_rsp).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_6 = k_vec.*(mu_TE.^2).*(omega_rsp).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                int_FCT_7 = k_vec.*(mu_TM.^2).*(omega_rsp).*(gamma./((w_mn{cb_num,vb_num}-w).^2+gamma^2));
                I_TE_FCT     = trapz(k_vec, int_FCT_1);
                I_TM_FCT     = trapz(k_vec, int_FCT_2);
                I_TE_N_FCT   = trapz(k_vec, int_FCT_3);
                I_TM_N_FCT   = trapz(k_vec, int_FCT_4);
                I_Rsp_FCT    = trapz(k_vec, int_FCT_5);
                I_Rsp_TE_FCT = trapz(k_vec, int_FCT_6);
                I_Rsp_TM_FCT = trapz(k_vec, int_FCT_7);
                                
                sum_TE_FCT = sum_TE_FCT + I_TE_FCT;
                sum_TM_FCT = sum_TM_FCT + I_TM_FCT;
                sum_TE_N_FCT = sum_TE_N_FCT + I_TE_N_FCT;
                sum_TM_N_FCT = sum_TM_N_FCT + I_TM_N_FCT;                
                sum_Rsp_FCT = sum_Rsp_FCT + I_Rsp_FCT;
                sum_Rsp_TE_FCT = sum_Rsp_TE_FCT + I_Rsp_TE_FCT;
                sum_Rsp_TM_FCT = sum_Rsp_TM_FCT + I_Rsp_TM_FCT;
                
%                 figure(h_FCT_Integrands);
%                 subplot(211);
%                 [AX,H1,H2] = plotyy(k_vec, real(int_FCT_1), k_vec, imag(int_FCT_1));
%                 set(get(AX(1), 'Ylabel'), 'String', '\Re(Int)');
%                 set(get(AX(2), 'Ylabel'), 'String', '\Im(Int)');
%                 %plot(k_vec, real(int_1), 'b', k_vec, imag(int_1), 'r');
%                 title(['FCT - cb=' num2str(cb_num) ', vb=' num2str(vb_num) ': \omega=' num2str(w, '%1.5e') ' s^{-1}']);
%                 subplot(212);
%                 [AX,H1,H2] = plotyy(k_vec, real(int_FCT_6), k_vec, imag(int_FCT_6));
%                 set(get(AX(1), 'Ylabel'), 'String', '\Re(Int)');
%                 set(get(AX(2), 'Ylabel'), 'String', '\Im(Int)');
%                 %plot(k_vec, real(int_3), 'b', k_vec, imag(int_3), 'r');
%                 xlabel('k [m^-^1]'); 
%                 drawnow;
                
                % HF
                [Q_k_TE, M_TE] = CalculateHFEnhancement(k_vec,omega_g, mu_TE, w_mn_tag{cb_num,vb_num}, w, gamma, Theta{cb_num,vb_num}, 'Exact');
                [Q_k_TM, M_TM] = CalculateHFEnhancement(k_vec,omega_g, mu_TM, w_mn_tag{cb_num,vb_num}, w, gamma, Theta{cb_num,vb_num}, 'Exact');
                M_Rsp = (2/3)*M_TE + (1/3)*M_TM;
                Q_k_Rsp = M_Rsp\ones(k_len,1);
                
                int_HF_1 = k_vec.*mu_TE.^2.*omega_g.*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(imag(Q_k_TE.').*(w_mn_tag{cb_num,vb_num}-w)+gamma.*real(Q_k_TE.'));
                int_HF_2 = k_vec.*mu_TM.^2.*omega_g.*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(imag(Q_k_TM.').*(w_mn_tag{cb_num,vb_num}-w)+gamma.*real(Q_k_TM.'));
                int_HF_3 = k_vec.*mu_TE.^2.*(omega_rsp).*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(imag(Q_k_TE.').*(w_mn_tag{cb_num,vb_num}-w)+gamma.*real(Q_k_TE.'));
                int_HF_4 = k_vec.*mu_TM.^2.*(omega_rsp).*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(imag(Q_k_TM.').*(w_mn_tag{cb_num,vb_num}-w)+gamma.*real(Q_k_TM.'));
                int_HF_5 = k_vec.*((2/3)*mu_TE.^2+(1/3)*mu_TM.^2).*(omega_rsp).*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(imag(Q_k_Rsp.').*(w_mn_tag{cb_num,vb_num}-w)+gamma.*real(Q_k_Rsp.'));
                int_HF_6 = k_vec.*mu_TE.^2.*(omega_g).*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(real(Q_k_TE.').*(w_mn_tag{cb_num,vb_num}-w)-gamma.*imag(Q_k_TE.'));
                int_HF_7 = k_vec.*mu_TM.^2.*(omega_g).*(-1./((w_mn_tag{cb_num,vb_num}-w).^2+gamma^2)).*(real(Q_k_TM.').*(w_mn_tag{cb_num,vb_num}-w)-gamma.*imag(Q_k_TM.'));
                I_TE_HF     = trapz(k_vec, int_HF_1);
                I_TM_HF     = trapz(k_vec, int_HF_2);
                I_TE_N_HF   = trapz(k_vec, int_HF_6);
                I_TM_N_HF   = trapz(k_vec, int_HF_7);
                I_Rsp_TE_HF = trapz(k_vec, int_HF_3);
                I_Rsp_TM_HF = trapz(k_vec, int_HF_4);
                %I_Rsp_TE_HF = trapz(k_vec, k_vec.*(mu_TE.^2.*Q_k_Rsp.').*(omega_rsp).*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)));
                %I_Rsp_TM_HF = trapz(k_vec, k_vec.*(mu_TM.^2.*Q_k_Rsp.').*(omega_rsp).*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)));
                I_Rsp_HF    = trapz(k_vec, int_HF_5);
                
                sum_TE_HF = sum_TE_HF + I_TE_HF;
                sum_TM_HF = sum_TM_HF + I_TM_HF;
                sum_TE_N_HF = sum_TE_N_HF + I_TE_N_HF;
                sum_TM_N_HF = sum_TM_N_HF + I_TM_N_HF;
                sum_Rsp_TE_HF = sum_Rsp_TE_HF + I_Rsp_TE_HF;
                sum_Rsp_TM_HF = sum_Rsp_TE_HF + I_Rsp_TM_HF;
                sum_Rsp_HF = sum_Rsp_HF + I_Rsp_HF;
                
%                 figure(h_HF_Integrands);
%                 subplot(211);
%                 [AX,H1,H2] = plotyy(k_vec, real(int_HF_1), k_vec, imag(int_HF_1));
%                 set(get(AX(1), 'Ylabel'), 'String', '\Re(Int)');
%                 set(get(AX(2), 'Ylabel'), 'String', '\Im(Int)');
%                 %plot(k_vec, real(int_1), 'b', k_vec, imag(int_1), 'r');
%                 title(['HF - cb=' num2str(cb_num) ', vb=' num2str(vb_num) ': E=' num2str(E_exc(ee)/Consts.e_0) ' eV']);
%                 subplot(212);
%                 [AX,H1,H2] = plotyy(k_vec, real(int_HF_3), k_vec, imag(int_HF_3));
%                 set(get(AX(1), 'Ylabel'), 'String', '\Re(Int)');
%                 set(get(AX(2), 'Ylabel'), 'String', '\Im(Int)');
%                 %plot(k_vec, real(int_3), 'b', k_vec, imag(int_3), 'r');
%                 xlabel('k [m^-^1]'); 
%                 drawnow;
            end
        end
     
        % FCT
        alpha_TE_FCT(ee)   = -C_g.*2*w.*sum_TE_FCT/100;       
        alpha_TM_FCT(ee)   = -C_g.*2*w.*sum_TM_FCT/100;
        delta_n_TE_FCT(ee) = -C_d_n.*sum_TE_N_FCT;
        delta_n_TM_FCT(ee) = -C_d_n.*sum_TM_N_FCT;
        Rsp_TE_FCT(ee)     = C_rsp.*w^3.*sum_Rsp_TE_FCT;
        Rsp_TM_FCT(ee)     = C_rsp.*w^3.*sum_Rsp_TM_FCT;
        Rsp_FCT(ee)        = (2/3)*Rsp_TE_FCT(ee) + (1/3)*Rsp_TM_FCT(ee); 
        
        % HF
        alpha_TE_HF(ee)   = -C_g*2*w*(sum_TE_HF)/100;
        alpha_TM_HF(ee)   = -C_g*2*w*(sum_TM_HF)/100;
        delta_n_TE_HF(ee) = -C_d_n*(sum_TE_HF);
        delta_n_TM_HF(ee) = -C_d_n*(sum_TM_HF);
        Rsp_TE_HF(ee)     = C_rsp.*w^3.*(sum_Rsp_TE_HF);
        Rsp_TM_HF(ee)     = C_rsp.*w^3.*(sum_Rsp_TM_HF);
        Rsp_HF(ee)        = (2/3)*Rsp_TE_HF(ee) + (1/3)*Rsp_TM_HF(ee); 

        D_E_CH_con = {D_E_CH_con, D_E_CH};
        D_E_SX_con = {D_E_SX_con, D_E_SX};
        
    end
    
    alpha_TE_FCT_vec{con_num} = alpha_TE_FCT;
    alpha_TM_FCT_vec{con_num} = alpha_TM_FCT;
    delta_n_TE_FCT_vec{con_num} = delta_n_TE_FCT;
    delta_n_TM_FCT_vec{con_num} = delta_n_TM_FCT;
    Rsp_FCT_vec{con_num} = Rsp_FCT;
    Rsp_TE_FCT_vec{con_num} = Rsp_TE_FCT;
    Rsp_TM_FCT_vec{con_num} = Rsp_TM_FCT;
    alpha_TE_HF_vec{con_num} = alpha_TE_HF;
    alpha_TM_HF_vec{con_num} = alpha_TM_HF;
    delta_n_TE_HF_vec{con_num} = delta_n_TE_HF;
    delta_n_TM_HF_vec{con_num} = delta_n_TM_HF;
    Rsp_HF_vec{con_num} = Rsp_HF;
    Rsp_TE_HF_vec{con_num} = Rsp_TE_HF;
    Rsp_TM_HF_vec{con_num} = Rsp_TM_HF;
    
    disp(['This stage took ', num2str(toc(t_start)/60), ' minutes']);
    disp(' ');
end

PlotAbsorptionSpEmission(E_exc, alpha_TE_FCT_vec, alpha_TM_FCT_vec, alpha_TE_HF_vec, alpha_TM_HF_vec, Rsp_FCT_vec, Rsp_HF_vec, delta_n_TE_FCT_vec,  delta_n_TM_FCT_vec,  delta_n_TE_HF_vec,  delta_n_TM_HF_vec, con_vec, T, E_f_vec, E_e_mat, E_h_mat, d_E, E_g);
