global Consts;

%% Part 1 - Subband energies
%close all; clc; clear all;
warning off;

%Init project
global project_path;
project_path = 'C:\Users\Yossi Mchaeli\Documents\Code\MATLAB\Thesis';
cd(project_path);
run('.\Common\AddPath.m');

Constants;

% well_width = 50   % [A]
% x = 0.2;
% Res1 = SingleQW_kp(well_width,x);

L_z = 50e-10;
x = 0.2;
Structure = { 'GaAlAs' , 100, x ;
              'GaAs', L_z*1e10, 0;
              'GaAlAs' , 100, x };

%Res2 = GeneralStructure_TransferMatrix(Structure);
%Res2 = BandStructureAndDispersion(Structure, 'Matrix');

close all;

%% Part 2 - Gain spectra

Con_vec = [2e11,6e11,2e12,3e12,4e12];    % [cm^-2]
T = 300;
init_E = 0.01*Consts.e_0;
Mb_s = sqrt(Consts.e_0*Res2.Materials{2}.E_p*Consts.m_0/6);

for (con_num = 1:length(Con_vec))
    N = Con_vec(con_num); P = N;
    
    options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
    E_fv = fzero(@(E_fv) TargetFermiIntegral(Res2.k_t_mat, Res2.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
    E_fc = fzero(@(E_fc) TargetFermiIntegral(Res2.k_t_mat, Res2.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);
    
    %E_exc = (Res2.Materials{2}.E_g-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (Res2.Materials{2}.E_g+0.5)*Consts.e_0;
    deltaE_vc_0 = (Res2.E_k_v{1}(1)/Consts.e_0 + Res2.E_k_c{1}(1)/Consts.e_0)
    E_exc = (deltaE_vc_0-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (deltaE_vc_0+0.5)*Consts.e_0;
    C = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*Res2.Materials{2}.n);
    gamma = 1e13;      % [s^-1]
    
    for (ii=1:length(E_exc))
        w = (E_exc(ii)./Consts.hbar);
        sum_TE = 0;
        sum_TM = 0;
        
        for (cb_num = 1:1)
            for (vb_num = 1:3)
                [ii, cb_num, vb_num]
                               
                mu_TE = sqrt(squeeze(Res2.rtsTE(vb_num,cb_num,:))).'.*Mb_s;
                mu_TM = sqrt(squeeze(Res2.rtsTM(vb_num,cb_num,:))).'.*Mb_s;
                f_h = FermiDirac(E_fv, Res2.E_k_v{vb_num}, T);
                f_e = FermiDirac(E_fc, Res2.E_k_c{cb_num}, T);
                w_mn = (Res2.E_k_v{vb_num} + Res2.E_k_c{cb_num})./Consts.hbar;
                I_TE = 2*trapz(Res2.k_t_mat, Res2.k_t_mat.*mu_TE.^2.*(f_e+f_h-1).*(gamma./((w_mn-w).^2+gamma^2)));
                I_TM = 2*trapz(Res2.k_t_mat, Res2.k_t_mat.*mu_TM.^2.*(f_e+f_h-1).*(gamma./((w_mn-w).^2+gamma^2)));
                
                sum_TE = sum_TE + I_TE;
                sum_TM = sum_TM + I_TM;
            end
        end
        
        G_TE(ii) = C.*w.*sum_TE;
        G_TM(ii) = C.*w.*sum_TM;
    end
    
    index_TE = find(abs(G_TE) == max(abs(G_TE(1:round(length(G_TE)/3)))));
    index_TM = find(abs(G_TM) == max(abs(G_TM(1:round(length(G_TM)/3)))));
    
    figure(1);
    subplot(211);
    plot(E_exc/Consts.e_0, G_TE/100); hold on;
    plot(deltaE_vc_0*ones(size(G_TE)), G_TE/100, 'g:');
    ylabel('G [cm^{-1}]');
    text(E_exc(index_TE)/Consts.e_0,G_TE(index_TE)/100, [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    subplot(212);
    plot(E_exc/Consts.e_0, G_TM/100); hold on;
    plot(deltaE_vc_0*ones(size(G_TE)), G_TE/100, 'g:');
    text(E_exc(index_TM)/Consts.e_0,G_TM(index_TM)/100, [num2str(N/1e12), 'x10^{12} cm^{-2}'], 'FontSize', 8);
    xlabel('E [eV]'); ylabel('G [cm^{-1}]');
end

%T_vec = [20:20:100, 100:50:300];
%Con = 2e11;

% for (T_num = 1:length(T_vec))
%     N = Con; P = N;
%     T = T_vec(T_num);
%     
%     options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
%     E_fv = fzero(@(E_fv) TargetFermiIntegral(Res2.k_t_mat, Res2.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
%     E_fc = fzero(@(E_fc) TargetFermiIntegral(Res2.k_t_mat, Res2.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);
%     
%     E_exc = (Res2.Materials{2}.E_g-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (Res2.Materials{2}.E_g+0.5)*Consts.e_0;
%     C = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*Res2.Materials{2}.n);
%     gamma = 1e13;      % [s^-1]
%     
%     for (ii=1:length(E_exc))
%         w = (E_exc(ii)./Consts.hbar);
%         sum_TE = 0;
%         sum_TM = 0;
%         
%         for (cb_num = 1:1)
%             for (vb_num = 1:3)
%                 [ii, cb_num, vb_num]
%                 
%                 %Mb_s = Consts.e_0*Res2.Materials{2}.E_p*Consts.m_0/6;
%                 mu_TE = sqrt(squeeze(Res2.rtsTE(vb_num,cb_num,:))).'.*Mb_s;
%                 mu_TM = sqrt(squeeze(Res2.rtsTM(vb_num,cb_num,:))).'.*Mb_s;
%                 f_h = FermiDirac(E_fv, Res2.E_k_v{vb_num}, T);
%                 f_e = FermiDirac(E_fc, Res2.E_k_c{cb_num}, T);
%                 w_mn = (Res2.E_k_v{vb_num} + Res2.E_k_c{cb_num})./Consts.hbar;
%                 I_TE = 2*trapz(Res2.k_t_mat, Res2.k_t_mat.*mu_TE.^2.*(f_e+f_h-1).*(gamma./((w_mn-w).^2+gamma^2)));
%                 I_TM = 2*trapz(Res2.k_t_mat, Res2.k_t_mat.*mu_TM.^2.*(f_e+f_h-1).*(gamma./((w_mn-w).^2+gamma^2)));
%                 
%                 sum_TE = sum_TE + I_TE;
%                 sum_TM = sum_TM + I_TM;
%             end
%         end
%         
%         G_TE(ii) = C.*w.*sum_TE*L_z;
%         G_TM(ii) = C.*w.*sum_TM*L_z;
%     end
%     
%     index_TE = find(abs(G_TE) == max(abs(G_TE(1:round(length(G_TE)/3)))));
%     index_TM = find(abs(G_TM) == max(abs(G_TM(1:round(length(G_TM)/3)))));
%     
%     figure(2);
%     subplot(211);
%     plot(E_exc/Consts.e_0, G_TE/100); hold on;
%     plot(Res2.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
%     ylabel('G [cm^{-1}]');
%     text(E_exc(index_TE)/Consts.e_0,G_TE(index_TE)/100, [num2str(T), 'K'], 'FontSize', 8);
%     subplot(212);
%     plot(E_exc/Consts.e_0, G_TM/100); hold on;
%     plot(Res2.Materials{2}.E_g*ones(size(G_TE)), G_TE/100, 'g:');
%     text(E_exc(index_TM)/Consts.e_0,G_TM(index_TM)/100, [num2str(T), 'K'], 'FontSize', 8);
%     xlabel('E [eV]'); ylabel('G [cm^{-1}]');
% end