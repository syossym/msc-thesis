function OpticalParams = CalculateOpticalParameters(QStruct, Cond, Valence, OpticalParams, params)
%
% This function calculates the optical susceptability of the structure
% using the Free-Carrier and Hartree-Fock model.
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'Cond' - structure containig the conduction band energies and
%                   wavefucntion.
%          'Valence' - structure containig the valence band energies and
%                      wavefucntion.
%          'OpticalParams' - structure containing pre-calculated optical
%                            parameters for the simulated structure.
%          'params' - simulaion parameters.
%
%   Output: -
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Simulation parameters

% Simulation constants
Mb = sqrt(Consts.e_0*QStruct.ActiveLayers{1}.E_p*Consts.m_0/6);
E_p = Consts.e_0*QStruct.ActiveLayers{1}.E_p;           % [J]
E_g = QStruct.ActiveLayers{1}.E_g;                      % [eV]
init_E = 0.001*Consts.e_0;
L_z = QStruct.ActiveLayers{1}.L*1e-10;                  % [m]
E_exc = params.E_exc;                                   % [eV]
T = params.T;   
gamma = params.gamma;

C_g = 1/(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*QStruct.ActiveLayers{1}.n);
C_d_n = 1/(Consts.hbar*Consts.eps_0*L_z*pi*QStruct.ActiveLayers{1}.n);
C_rsp = QStruct.ActiveLayers{1}.n/(pi^3*Consts.hbar^2*Consts.c^3*L_z*Consts.eps_0);

% Simulation grids 
sim_index = params.k_indices; %length(R.k_t_mat);
dk = QStruct.k_t_vec(2) - QStruct.k_t_vec(1);
k_vec = interp(QStruct.k_t_vec(sim_index), params.k_scale);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.01:2*pi];

%% Variable pre-processing 

% Reshaping the data
for (cb_num = 1:params.num_cond_subbands)
    for (vb_num = 1:params.num_valence_subbands)
        E_e_mat{cb_num} = interp1(QStruct.k_t_vec(sim_index), Cond.E_k(cb_num,sim_index), k_vec, 'pchip');
        E_h_mat{vb_num} = interp1(QStruct.k_t_vec(sim_index), Valence.E_k(vb_num,sim_index), k_vec, 'pchip');
        mu_TE_mat{cb_num,vb_num} = interp1(QStruct.k_t_vec(sim_index), squeeze(OpticalParams.rts_TE(vb_num,cb_num,sim_index)), k_vec, 'pchip');
        mu_TM_mat{cb_num,vb_num} = interp1(QStruct.k_t_vec(sim_index), squeeze(OpticalParams.rts_TE(vb_num,cb_num,sim_index)), k_vec, 'pchip');
    end
end

% Form factor calculation 
disp('-- Calculating form factors');
for (cb_num = 1:params.num_cond_subbands)
    for (vb_num = 1:params.num_valence_subbands)
        disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);

        Psi_e = abs(Cond.Wf_0_Q(cb_num,:)).';
        Psi_hh = abs(Valence.Wf_0_Q{vb_num}(:,2));
        Psi_lh = abs(Valence.Wf_0_Q{vb_num}(:,3));

        Xi_e = repmat(Psi_e.^2, 1, length(QStruct.ActiveLayers{1}.z_grid));
        Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(QStruct.ActiveLayers{1}.z_grid), 1);
        Xi_ee = Xi_e.*Xi_e.';
        Xi_hh = Xi_h.'.*Xi_h;
        Xi_eh = Xi_e.*Xi_h;

        Gq_eh{cb_num,vb_num} = G_q(Xi_eh,k_vec,QStruct.ActiveLayers{1}.z_grid);
        Gq_ee{cb_num,vb_num} = G_q(Xi_ee,k_vec,QStruct.ActiveLayers{1}.z_grid);
        Gq_hh{cb_num,vb_num} = G_q(Xi_hh,k_vec,QStruct.ActiveLayers{1}.z_grid);
    end
end


%% Simulation parameters calculation

t_start = tic;

% Carrier concentrations
P = params.P;                      % [cm^-2]
N = params.N_DEG + P;

disp(['-- Concentration - N=' num2str(N, '%1.0e') ' cm^-2, P=' num2str(P, '%1.0e') ' cm^-2']);

%E_fc = 3.6e-11*N*1e-3*Consts.e_0 + R.E_k_c{1}(1);

options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000);
E_fv = fzero(@(E_fv) TargetFermiIntegral(k_vec, E_h_mat, E_fv, P*1e4, T, L_z, 'h'), init_E, options);
E_fc = fzero(@(E_fc) TargetFermiIntegral(k_vec, E_e_mat, E_fc, N*1e4, T, L_z, 'c'), init_E, options);

E_f_vec = E_fc;

%PlotQuasiFermiLevel(h_Band_Diagrams, E_fc, E_fv, k_vec, T, E_e_mat, E_h_mat, N, P);
%PlotFermiLevel(h_Band_Diagrams, E_fc,  QStruct.k_t_vec, T, R.E_k_c, R.E_k_v, N);

D_E_CH_con = {};
D_E_SX_con = {};
w_tag_con  = {};

disp('Calculating screened dielectric function');
[eps_q, eps_q_fit] = CalculateScreening(k_vec, theta, E_fv, E_fc, T, E_e_mat(1:params.num_cond_subbands),E_h_mat(1:params.num_valence_subbands), QStruct.ActiveLayers{1}.eps_r, Gq_ee, Gq_hh);

for (cb_num = 1:params.num_cond_subbands)
    for (vb_num = 1:params.num_valence_subbands)
        disp([' - cb=' num2str(cb_num) ', vb=' num2str(vb_num)]);

        E_e = E_e_mat{cb_num};
        E_h = E_h_mat{vb_num};

        %f_h4 = 1 - exp(-(E_h-E_fc)/(Consts.k_B*T));
        f_h = FermiDirac(E_fv, E_h, T);
        f_e = FermiDirac(E_fc, E_e, T);

        disp('Calculating Theta matrix and Screened-exchange shift energy');
        [Theta{cb_num,vb_num}, D_E_SX{cb_num,vb_num}] = CalculateThetaESX(Gq_ee{cb_num,vb_num},Gq_hh{cb_num,vb_num},Gq_eh{cb_num,vb_num},theta,eps_q_fit,k_vec,f_e,f_h,QStruct.ActiveLayers{1}.eps_r);

        disp('Calculating D_E_CH');
        D_E_CH{cb_num,vb_num} = CalculateECH(Gq_hh{cb_num,vb_num},eps_q_fit,k_vec,QStruct.ActiveLayers{1}.eps_r);

        disp('Calculating bandgaps');
        w_mn_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num})./Consts.hbar;
        %w_mn_tag{cb_num,vb_num} = (E_h + E_e + D_E_CH{cb_num,vb_num}.*ones(size(k_vec)))./Consts.hbar;
        w_mn{cb_num,vb_num} = (E_h + E_e)./Consts.hbar;

        delta_renorm{cb_num,vb_num} = (D_E_CH{cb_num,vb_num}.*ones(size(k_vec)) + D_E_SX{cb_num,vb_num});
        %delta_renorm{cb_num,vb_num} = (D_E_CH{cb_num,vb_num}.*ones(size(k_vec)));
        d_E(1,cb_num,vb_num) = delta_renorm{cb_num,vb_num}(1);
    end
end

%% Spectrum calculation

h_wait = waitbar(0,['Calculating Spectrum, N_{DEG}=' num2str(params.N_DEG,'%1.0e') 'cm^{-2}']);

%h_FCT_Integrands = figure;
%h_HF_Integrands = figure;
for (ee = 1:length(E_exc))
    waitbar(ee/length(E_exc),h_wait);

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

    for (cb_num = 1:params.num_cond_subbands)
        for (vb_num = 1:params.num_valence_subbands)

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

            % HF
            [Q_k_TE, M_TE] = CalculateHFEnhancement(k_vec,omega_g, mu_TE, w_mn_tag{cb_num,vb_num}, w, gamma, Theta{cb_num,vb_num}, 'Exact');
            [Q_k_TM, M_TM] = CalculateHFEnhancement(k_vec,omega_g, mu_TM, w_mn_tag{cb_num,vb_num}, w, gamma, Theta{cb_num,vb_num}, 'Exact');
            M_Rsp = (2/3)*M_TE + (1/3)*M_TM;
            Q_k_Rsp = M_Rsp\ones(k_len,1);

            int_HF_1 = k_vec.*mu_TE.^2.*omega_g.*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TE.';
            int_HF_2 = k_vec.*mu_TM.^2.*omega_g.*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_TM.';
            int_HF_3 = k_vec.*(mu_TE.^2.*Q_k_TE.').*(omega_rsp).*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma));
            int_HF_4 = k_vec.*(mu_TM.^2.*Q_k_TM.').*(omega_rsp).*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma));
            int_HF_5 = k_vec.*((2/3)*mu_TE.^2+(1/3)*mu_TM.^2).*(omega_rsp).*(1i./(1i*(w_mn_tag{cb_num,vb_num}-w)+gamma)).*Q_k_Rsp.';
            I_TE_HF     = trapz(k_vec, int_HF_1);
            I_TM_HF     = trapz(k_vec, int_HF_2);
            I_TE_N_HF   = trapz(k_vec, int_HF_1);
            I_TM_N_HF   = trapz(k_vec, int_HF_2);
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
    Xi_TE_FCT(ee)      = (2/QStruct.ActiveLayers{1}.n)*delta_n_TE_FCT(ee)+1i*(Consts.c/w/QStruct.ActiveLayers{1}.n)*(alpha_TE_FCT(ee)*100);
    Xi_TM_FCT(ee)      = (2/QStruct.ActiveLayers{1}.n)*delta_n_TM_FCT(ee)+1i*(Consts.c/w/QStruct.ActiveLayers{1}.n)*(alpha_TM_FCT(ee)*100);

    % HF
    alpha_TE_HF(ee)   = -C_g*2*w*imag(sum_TE_HF)/100;
    alpha_TM_HF(ee)   = -C_g*2*w*imag(sum_TM_HF)/100;
    delta_n_TE_HF(ee) = -C_d_n*real(sum_TE_HF);
    delta_n_TM_HF(ee) = -C_d_n*real(sum_TM_HF);
    Rsp_TE_HF(ee)     = C_rsp.*w^3.*imag(sum_Rsp_TE_HF);
    Rsp_TM_HF(ee)     = C_rsp.*w^3.*imag(sum_Rsp_TM_HF);
    Rsp_HF(ee)        = (2/3)*Rsp_TE_HF(ee) + (1/3)*Rsp_TM_HF(ee);
    Xi_TE_HF(ee)      = (2/QStruct.ActiveLayers{1}.n)*delta_n_TE_HF(ee)+1i*(Consts.c/w/QStruct.ActiveLayers{1}.n)*(alpha_TE_HF(ee)*100);
    Xi_TM_HF(ee)      = (2/QStruct.ActiveLayers{1}.n)*delta_n_TM_HF(ee)+1i*(Consts.c/w/QStruct.ActiveLayers{1}.n)*(alpha_TM_HF(ee)*100);

    D_E_CH_con = {D_E_CH_con, D_E_CH};
    D_E_SX_con = {D_E_SX_con, D_E_SX};

end
close(h_wait);

disp(['This stage took ', num2str(toc(t_start)/60), ' minutes']);
disp(' ');

%% Saving the results

OpticalParams.alpha_TE_FCT_vec = alpha_TE_FCT;
OpticalParams.alpha_TM_FCT_vec = alpha_TM_FCT;
OpticalParams.delta_n_TE_FCT_vec = delta_n_TE_FCT;
OpticalParams.delta_n_TM_FCT_vec = delta_n_TM_FCT;
OpticalParams.Rsp_FCT_vec = Rsp_FCT;
OpticalParams.Rsp_TE_FCT_vec = Rsp_TE_FCT;
OpticalParams.Rsp_TM_FCT_vec = Rsp_TM_FCT;
OpticalParams.Xi_TE_FCT_vec = Xi_TE_FCT;
OpticalParams.Xi_TM_FCT_vec = Xi_TM_FCT;

OpticalParams.alpha_TE_HF_vec = alpha_TE_HF;
OpticalParams.alpha_TM_HF_vec = alpha_TM_HF;
OpticalParams.delta_n_TE_HF_vec = delta_n_TE_HF;
OpticalParams.delta_n_TM_HF_vec = delta_n_TM_HF;
OpticalParams.Rsp_HF_vec = Rsp_HF;
OpticalParams.Rsp_TE_HF_vec = Rsp_TE_HF;
OpticalParams.Rsp_TM_HF_vec = Rsp_TM_HF;
OpticalParams.Xi_TE_HF_vec = Xi_TE_HF;
OpticalParams.Xi_TM_HF_vec = Xi_TM_HF;
