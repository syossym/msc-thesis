
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    This script calculates for a single GaAs/AlGaAs QW:
%    1. Conduction band dispersion + envelope function using the 
%       the shooting method.
%    2. Valence band dispersion + envelope functions using diagonalized 
%       4X4 k.p or 8X8 k.p method. For the 8X8 case the conduction bands
%       data used is obtained from the k.p method and not using the
%       shhoting method.
%    3. Optical transition strength (optical momentum matrix ratio between
%       the QW and the bulk semicondictor) between the conduction band and
%       the valence bands.
%    4. DOS for the valence bands calculated as 
%            dos_E = (1/pi)*(dk/dE)
%    5. Gain spectrum with zero and finite linewidths.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_cell = -1;    % marks the number of the cell to run. 
                  % Choose '-1' for the execution of the entire script.

if (run_cell == -1)
    clear all; close all; clc;
    run_cell = -1;
end
                  
%% 1. Constants (all MKS, except energy which is in eV)

if (run_cell == 1 || run_cell == -1)
   
    hbar = 1.055e-34;
    q = 1.602e-19;
    a = 3e-10;
    m_0 = 9.110e-31;
    c = 2.99792461e+8;     % speed of light
    e_0 = 8.854187827e-12; % permittivity of free space
    k_B = 1.380662e-23;    % Boltzmann constant
    Ep = 25.5;

end

%% 2. Plotting definitions

if (run_cell == 2 || run_cell == -1)

    h_struct = figure('Name','Structure');
    h_env_func = figure('Name','Envelope Functions');
    h_disp = figure('Name', 'Dispersion');
    h_rts = figure('Name', 'Transition Matrix');
    h_dos = figure('Name', 'DOS');
    h_gain = figure('Name', 'Gain');
    h_abs = figure('Name', 'Absorption');

end

%% 3. Luttinger-Kohn parameters

if (run_cell == 3 || run_cell == -1)

    g1 = 6.85; g2 = 2.1; g3 = 2.9;    %GaAs

    Gw1 = (hbar^2)*g1/(2*m_0);
    Gw2 = (hbar^2)*g2/(2*m_0);
    Gw3 = (hbar^2)*g3/(2*m_0);
    Cw = Gw1+Gw2;
    Dw = -2*Gw2+Gw1;
    Aw = Gw1-Gw2;
    Bw = 2*Gw2+Gw1;
    Ew = Gw2+Gw3;

    g1 = 3.45; g2 = 0.68; g3 = 1.29;  %AlAs

    a1 = (hbar^2)*g1/(2*m_0);
    Gb1 = (.9*Gw1)+(.1*a1);
    a2 = (hbar^2)*g2/(2*m_0);
    Gb2 = (.9*Gw2)+(.1*a2);
    a3 = (hbar^2)*g3/(2*m_0);
    Gb3 = (.9*Gw3)+(.1*a3);
    Cb = Gb1+Gb2;
    Db = -2*Gb2+Gb1;
    Ab = Gb1-Gb2;
    Bb = 2*Gb2+Gb1;
    Eb = Gb2+Gb3;

    Ev=0; Evb=((0.9*0)+(0.1*0.75))*q;
    Eg = 1.426;   % eV

end

%% 5. Structure parameters

if (run_cell == 5 || run_cell == -1)

    Nw = 70;
    Nb = Nw;
    Np = Nb+Nw+Nb;
    W = (Nw-1)*a*1e9;
    Z=zeros(Np);

    z = linspace(0,Np*a,Np);
    dz = a;
    dzz = dz^2;

    A  = [Ab*ones(1,Nb),Aw*ones(1,Nw),Ab*ones(1,Nb)];
    B  = [Bb*ones(1,Nb),Bw*ones(1,Nw),Bb*ones(1,Nb)];
    C  = [Cb*ones(1,Nb),Cw*ones(1,Nw),Cb*ones(1,Nb)];
    D  = [Db*ones(1,Nb),Dw*ones(1,Nw),Db*ones(1,Nb)];
    E  = [Eb*ones(1,Nb),Ew*ones(1,Nw),Eb*ones(1,Nb)];
    G2 = [Gb2*ones(1,Nb),Gw2*ones(1,Nw),Gb2*ones(1,Nb)];
    G3 = [Gb3*ones(1,Nb),Gw3*ones(1,Nw),Gb3*ones(1,Nb)];
    V  = [Evb*ones(1,Nb),zeros(1,Nw), Evb*ones(1,Nb)];

end

%% 6. Conduction band (using the shooting technique)

if (run_cell == 6 || run_cell == -1)

    hs_structure = [Nb*a, Nb, 0.1 ; (W/1e9), Nw, 0 ; Nb*a, Nb, 0.1];
    [E_c,V_cb,V_vb,m_eff,con_profile] = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
    wf_c  = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
    wf_c  = wf_c./sqrt(trapz(z, wf_c(1,:).^2));
    Evb = max(V_vb(:,2))./q;
    V_vb = -V_vb(:,2)./q;

    % Plot the structure -----------------------------------------------
    figure(h_struct);
    plot(z/1e-10, Eg*ones(1,Np) + V_cb(:,2)'./q, 'g', z/1e-10, V_vb, 'g', 'linewidth', 2);

end

%% 7. Velence band (4X4 k.p)

if (run_cell == 7 || run_cell == -1)

    for nk=1:100
        k(nk)=(nk-1)/500; % in A^-1
        l=0;m=1;lm=sqrt((l^2)+(m^2));
        kx=(l/lm)*k(nk)*1e10;ky=(m/lm)*k(nk)*1e10;
        k2=(kx^2)+(ky^2);
        k_t = sqrt(k2);
        k_t_vec(nk) = k_t;

        % Building the 2X2 diagonalized k.p Hamiltonian matrix ---------
        H = zeros(2*Np);
        for (zz = 1:Np)
            if (zz==1)
                H_diag  = [C(zz)*k_t^2+(D(zz+1)+2*D(zz))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                 ;
                                    sqrt(3)*G2(zz)*k_t^2               , A(zz)*k_t^2+(B(zz+1)+2*B(zz))/(2*dzz)+V(zz)   ];

                H_off_p = [       -(D(zz+1)+D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ;
                            sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)        ];

                H_off_m = [      -(D(zz))/(2*dzz)        ,      sqrt(3)*k_t*(G3(zz))/(2*dz)       ;
                           -sqrt(3)*k_t*(G3(zz))/(2*dz)  ,            -(B(zz))/(2*dzz)            ];

                H(1:2,1:2) = H_diag;
                H(1:2,3:4) = H_off_p;
                H(3:4,1:2) = H_off_m;
            elseif (zz==Np)
                H_diag  = [C(zz)*k_t^2+(2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                       ;
                                    sqrt(3)*G2(zz)*k_t^2               , A(zz)*k_t^2+(2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];

                H_off_p = [       -(D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz))/(2*dz)  ;
                            sqrt(3)*k_t*(G3(zz))/(2*dz)   ,            -(B(zz))/(2*dzz)        ];

                H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)        ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                           -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz) ,            -(B(zz)+B(zz-1))/(2*dzz)            ];

                H(2*Np-1:2*Np,2*Np-1:2*Np) = H_diag;
            else
                H_diag  = [C(zz)*k_t^2+(D(zz+1)+2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                       ;
                                    sqrt(3)*G2(zz)*k_t^2                       , A(zz)*k_t^2+(B(zz+1)+2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];

                H_off_p = [       -(D(zz+1)+D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ;
                            sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)        ];

                H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)         ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                           -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ,            -(B(zz)+B(zz-1))/(2*dzz)            ];

                H(2*zz-1:2*zz,2*zz-1:2*zz)   = H_diag;
                H(2*zz-1:2*zz,2*zz+1:2*zz+2) = H_off_p;
                H(2*zz+1:2*zz+2,2*zz-1:2*zz) = H_off_m;
            end
        end

        % Calculate the eigenvectors of the Hamiltonian ----------------
        [nk sum(sum(abs(H-H')))];
        [F,En] = eig(H);
        En = diag(En)./q;
        [En,I] = sort(real(En));
        En = -En;
        F1{nk} = F(1:2:end,I);
        F2{nk} = F(2:2:end,I);

        E_v{nk} = En(1:4);

        E1(nk)=En(1); E2(nk)=En(2); E3(nk)=En(3); E4(nk)=En(4);
        % E5(nk)=En(5);E6(nk)=En(6);E7(nk)=En(7);E8(nk)=En(8);
        % E9(nk)=En(9);E10(nk)=En(10);E11(nk)=En(11);E12(nk)=En(12);

        for (ii=1:length(E_v{1}))
            F_U_norm = sqrt(trapz(z, F1{nk}(:,ii).^2+F2{nk}(:,ii).^2));
            Fhh{nk}(:,ii) = F1{nk}(:,ii)./F_U_norm;
            Flh{nk}(:,ii) = F2{nk}(:,ii)./F_U_norm;
        end

        % Overlap integrals --------------------------------------------
        for (cb_index=1:length(E_c))
            for (vb_index=1:length(E_v{1}))
                overlap1(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*Fhh{nk}(:,vb_index)).^2;
                overlap2(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*Flh{nk}(:,vb_index)).^2;

                rtsTE(nk,cb_index,vb_index) = (3/4)*overlap1(nk, cb_index,vb_index) + ...
                    (3/12)*overlap2(nk, cb_index,vb_index);
                rtsTE_k(vb_index,cb_index,nk) = rtsTE(nk,cb_index,vb_index);
                rtsTM(nk,cb_index,vb_index) = overlap2(nk, cb_index,vb_index);
                rtsTM_k(vb_index,cb_index,nk) = rtsTM(nk,cb_index,vb_index);
            end
        end

        % Plot the envelope functions ----------------------------------
        figure(h_env_func);
        clf;
        colors = [1 0 0; 1 1 1; 0 1 0; 0 1 1];
        for(ii=1:length(E_v{1}))
            subplot(211); grid on; box on;
            hold on; h_F1(ii) = plot(z/1e-10, Fhh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
            l_F1{ii} = sprintf('F_h_h - h%d',ii);  ylabel('F_h_h');
            ylabel('F_h_h');
            subplot(212); grid on; box on;
            hold on; h_F2(ii) = plot(z/1e-10, Flh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
            l_F2{ii} = sprintf('F_l_h - h%d',ii); ylabel('F_l_h');
            ylabel('F_l_h');
        end
        subplot(211);
        hold on; plot(z/1e-10, wf_c, 'r--', 'linewidth', 2); hold off;
        legend(h_F1,l_F1); xlabel('z (A)');
        title(['Envelope function amplitudes, k_t=',num2str(k_t/1e10),' A^-^1']);
        subplot(212);
        hold on; plot(z/1e-10, wf_c, 'r--', 'linewidth', 2); hold off;
        legend(h_F2,l_F2); xlabel('z (A)');

    end

end

%% 9. DOS

if (run_cell == 9 || run_cell == -1)

    E_grid_step = 0.1e-3;   % eV
    k_t_step = k_t_vec(2)-k_t_vec(1);

    E_v_k = zeros(length(E_v{1}), length(k_t_vec));
    for (vb_index=1:length(E_v{1}))
        for (nk = 1:length(k_t_vec))
            E_v_k(vb_index, nk) = abs(E_v{nk}(vb_index));   % eV
        end
    end

    E_min = min(abs(V_vb));     % eV
    E_max = max(abs(V_vb));     % eV
    E_grid = [E_min:E_grid_step:E_max];    % energy grid for the DOS (eV)
    E_grid_size = size(E_grid,2);
    Dos_E = zeros(length(E_v{1}), length(E_grid));
    dE_k = diff(E_v_k,1,2);     % eV
    dk_t = diff(k_t_vec);       % m
    k_t_mid = k_t_vec(2:length(k_t_vec)) - k_t_step/2;
    E_v_k_node_index = (E_v_k-E_min)./E_grid_step + 1;

    for (vb_index = 1:length(E_v{1}))
        Dos_k(vb_index,:) = abs(k_t_mid.*dk_t./(pi*dE_k(vb_index,:)));  % 1/(eV m^2)
        for (step = 1:size(dE_k,2))
            E_1 = min(E_v_k(vb_index,step), E_v_k(vb_index,step+1));
            E_2 = max(E_v_k(vb_index,step), E_v_k(vb_index,step+1));
            node_index_1 = floor((E_1-E_min)/E_grid_step+1);
            node_index_2 = ceil((E_2-E_min)/E_grid_step+1);
            if (node_index_2<E_grid_size)
                for (node_index = node_index_1:(node_index_2-1))
                    Dos_E(vb_index, node_index) = Dos_E(vb_index, node_index) + Dos_k(vb_index, step);
                end
                node_fraction_1 = (E_1-E_grid(node_index_1))/E_grid_step;
                Dos_E(vb_index,node_index_1) = Dos_E(vb_index,node_index_1) - node_fraction_1*Dos_k(vb_index,step);
                node_fraction_2 = (E_grid(node_index_2)-E_2)/E_grid_step;
                Dos_E(vb_index,node_index_2-1) = Dos_E(vb_index,node_index_2-1) - node_fraction_2*Dos_k(vb_index,step);
            end
        end
    end

    diff_Dos_k = diff(Dos_k,1,2);
    Dos_k = [Dos_k(:,1), Dos_k(:,1:length(k_t_vec)-2)+diff_Dos_k, Dos_k(:,length(k_t_vec)-1)];

end

%% 10. Gain 

if (run_cell == 10 || run_cell == -1)
    
    Constants;
    
    % Parameters -------------------------------------------------------
    N_v = 2e1;           % valence band population (1e10 cm^-2)
    N_c = 2e1;             % conduction band population (1e10 cm^-2)
    Gamma = 0.02*q;      % Lorentzian linewidth
    T = 2;             % temperature (K)
    vb_index = 1;
    cb_index = 1;
    
    E_grid_step = 0.01*1e-3*q;      % J
    E_pump = [1.4*q:E_grid_step:2*q];    % grid for the pump energy
    % in excess of band gap energy (J)
    
    alpha_E_TE_tot = zeros(length(E_pump), 1);
    alpha_E_TM_tot = zeros(length(E_pump), 1);
    
    for(vb_index = 1:1)
        for (cb_index = 1:1)
            [ vb_index , cb_index ]
            
            alpha_E_TE = zeros(length(E_pump), 1);
            alpha_E_TM = zeros(length(E_pump), 1);
            
            C_0 = (pi*q^2*hbar)/(3.3*c*e_0*m_0^2);
            Mb_s = q*Ep*m_0/6;
            M_TE = Mb_s*squeeze(abs(rtsTE(:,cb_index,vb_index)));
            M_TM = Mb_s*squeeze(abs(rtsTM(:,cb_index,vb_index)));
            E_v_curr = -squeeze(E_v_k(vb_index,:)).*q;
            E_c_curr = 1.52*q + E_c(cb_index)*q + hbar^2*k_t_vec.^2/(2*0.0665*m_0);
            del_E_curr = E_c_curr - E_v_curr;
            
            dos_v = Dos_E(vb_index, Dos_E(vb_index,:) ~= 0); dos_v = dos_v(1)/q;
            Ef_v = -QuasiFermiLevels(max(abs(V_vb)), -E_v_curr(1), dos_v, N_v*1e10*1e4, T);
            dos_c = 1e-4*0.0665*m_0/(pi*hbar^2);                                         % conduction band density of states [cm^-2]
            Ef_c = 1.52*q + E_c(cb_index)*q + k_B*T*log(exp(1e10*N_c/(dos_c*k_B*T))-1);  % conduction subband QF level, relative to the band edge
            f_v = FermiDirac(Ef_v, E_v_curr, T);
            f_c = FermiDirac(Ef_c, E_c_curr, T);
            del_f = f_c - f_v;
            
            for (aa = 1:length(E_pump))
                
                E_curr = E_pump(aa);
                B = (Gamma/(2*pi))./((del_E_curr-E_curr).^2+0.25*Gamma^2);
                alpha_E_TE(aa) = 4*C_0*(1/E_curr)*trapz(k_t_vec, (k_t_vec/(2*pi)).*M_TE'.*del_f.*B )/(pi*W*1e-9);
                alpha_E_TM(aa) = 4*C_0*(1/E_curr)*trapz(k_t_vec, (k_t_vec/(2*pi)).*M_TM'.*del_f.*B )/(pi*W*1e-9);
                
            end
            
            alpha_E_TE_tot = alpha_E_TE_tot + alpha_E_TE;
            alpha_E_TM_tot = alpha_E_TM_tot + alpha_E_TM;
        end
    end
    
    figure(1);
    plot(E_pump/(1e-3*q), 1e-2*alpha_E_TE_tot, 'b', E_pump/(1e-3*q), 1e-2*alpha_E_TM_tot, 'g');
    legend('TE', 'TM'); xlabel('E_{pump} [meV]'); ylabel('\alpha(E_{pump}) [cm^-^1]');
    grid on; box on;
        
        %     E_grid_step = 0.1*1e-3*q;      % J
        %
        %     E_pump = [0:E_grid_step:500*1e-3*q];    % grid for the pump energy
        %                                             % in excess of band gap energy (J)
    %     gain_E = zeros(length(E_pump), 1);
    %
%     % Gain calculation for all valence & conduction combinations -------
%     for (vb_index = 1:length(E_v{1}))
%         for (cb_index=1:length(E_c))
%             
%             [vb_index, cb_index]
%             
%             rts = squeeze(rtsTE_k(vb_index,cb_index,:));
%             Dos_k_spec = squeeze(Dos_k(vb_index,:))/q;    % 1/(J m^2) 
%             E_n = squeeze(E_v_k(vb_index,:)).*q;
% 
%             gain_k = (q^2)*hbar/(e_0*c*m_0^2)*rts'./(1./Dos_k_spec + pi*hbar^2/(0.0665*m_0));   % g(k)
%             gain_k = gain_k./(1.52*q + E_n + E_c(cb_index)*q + hbar^2*k_t_vec.^2/(2*0.0665*m_0));
%             gain_k = gain_k*28.8*q*m_0/2;
% 
%             I = find(Dos_E(vb_index, :) ~= 0);
%             dE_f_v = QuasiFermiLevels((V_vb)*q, (E_v{1}(vb_index))*q, Dos_E(vb_index, I(1))/q, 1e4*N_v*1e10, T);
%             E_f_v = abs(E_v{1}(vb_index))*q + dE_f_v;
%             
%             % QuasiFermiLevels(vb_index, N_v, T);         % valence subband quasi-Fermi level, relative to the valence band edge
% 
%             Dos_c = 1e-4*m_eff(Np/2,2)/(pi*hbar^2);               % density of states in the conduction band (cm^-2 J^-1)
%             dE_f_c = k_B*T*log(exp(1e10*N_c/(Dos_c*k_B*T))-1);    % conduction subband quasi-Fermi level, relative to the conduction subband edge
%             E_f_c = E_c(cb_index)*q + dE_f_c;                     % conduction band quasi-Fermi level, relative to the conduction band edge
% 
%             gain_k = gain_k.*(FermiDirac(E_f_v, E_v_k(vb_index,:)*q, T) + ...
%                              FermiDirac(E_f_c, E_c(cb_index)*q+hbar^2*k_t_vec.^2/(2*0.0665*m_0), T) - ...
%                               1);
%             E_diff = E_v_k(vb_index,:)*q + E_c(cb_index)*q + hbar^2*k_t_vec.^2/(2*0.0665*m_0);  % energy difference between the conduction and the valence bands vs. k
%             dE_diff = diff(E_diff);
% 
%             % Fitting Gain(k) to Gain(E)
%             for (step=1:size(dE_diff,2))
%                 E_1 = min(E_diff(step), E_diff(step+1));
%                 E_2 = max(E_diff(step), E_diff(step+1));
% 
%                 node_index_1 = floor(E_1/E_grid_step+1);
%                 node_index_2 = ceil(E_2/E_grid_step+1);
%                 if (node_index_2 < length(E_pump))
%                     for (node_index=node_index_1:(node_index_2-1))
%                         gain_E(node_index) = gain_E(node_index) + gain_k(step);
%                     end
%                     node_fraction_1 = (E_1-E_pump(node_index_1))/E_grid_step;
%                     gain_E(node_index_1) = gain_E(node_index_1) - node_fraction_1*gain_k(step);
%                     node_fraction_2 = (E_pump(node_index_2)-E_2)/E_grid_step;
%                     gain_E(node_index_2-1) = gain_E(node_index_2-1) - node_fraction_2*gain_k(step);
%                 end
%             end
% 
%             % Adding linewidth
%             if (Gamma ~= 0)
%                 factor = 1/0.844; % this factor maekes up for the loss due to
%                 % limiting line shape width to 2*Gamma
% 
%                 gain_E_lw = zeros(length(E_pump), 1);
%                 width = 2*Gamma;
%                 dE = E_grid_step;
%                 E_width = floor(width/dE);
% 
%                 for (ii=1:length(E_pump))
%                     for (jj=ii-E_width:ii+E_width)
%                         if (jj>0 && jj<length(E_pump))
%                             fract = dE*(Gamma/(2*pi))/((E_pump(ii)-E_pump(jj))^2+(Gamma/2)^2);
%                             gain_E_lw(jj) = gain_E_lw(jj) + fract + gain_E(ii);
%                         end
%                     end
%                 end
% 
%                 gain_E = gain_E_lw*factor;
%             end
% 
%             Gain_E{vb_index, cb_index, :} = gain_E;
%             Gain_k{vb_index, cb_index, :} = gain_k;
%         end
%     end

end

%% 11. Absorption 

if (run_cell == 11) %|| run_cell == -1)

    sums_TE = [];
    sums_TM = [];
    
    Gamma = 2e-3*q;           % J
    w_grid = k_t_vec*c;
    w_grid(1) = w_grid(2);
    E_grid = w_grid*hbar;     % J
    
    rts_TE = abs(squeeze(rtsTE_k(:,1,:)));
    rts_TM = abs(squeeze(rtsTM_k(:,1,:)));
    
    for (vb_index = 1:length(E_v{1}))   
        del_E = (E_c(1) + Eg + E_v{1}(vb_index))*q;    % J
        line_width = (Gamma/(2*pi))./((del_E*ones(size(E_grid))-E_grid).^2+ones(size(E_grid))*(Gamma/2)^2);
        sums_TE = [sums_TE; rts_TE(vb_index,:).*line_width];
        sums_TM = [sums_TM; rts_TM(vb_index,:).*line_width];
    end

    alpha_TE = (1./w_grid) .* sum(sums_TE);
    alpha_TM = (1./w_grid) .* sum(sums_TM);

end

%% 12. Plotting  

if (run_cell == 12 || run_cell == -1)

    % k_t=0 bands and envelope functions -------------------------------
    figure(h_struct); grid on;
    for (cb_index=1:length(E_c))
        hold on;
        plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, Eg + E_c(cb_index)*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), 'r', 'linewidth', 2 );
        text(z(Nb+Nw+10)/1e-10, Eg + E_c(cb_index), sprintf('e%d',cb_index));
        plot(z/1e-10, (Eg + E_c(cb_index))*ones(1,Np) + wf_c(cb_index,:)./1e5, 'r');
    end
    for (vb_index=1:length(E_v{1}))
        hold on;
        plot(z(V_vb==max(V_vb))/1e-10, E_v{1}(vb_index)*ones(1,length(z(V_vb==max(V_vb)))), 'b', 'linewidth', 2);
        text(z(Nb+Nw+10)/1e-10, E_v{1}(vb_index), sprintf('h%d',vb_index));
        plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np)*q + Fhh{1}(:,vb_index).'./1e5, ':b');
        plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np)*q + Flh{1}(:,vb_index).'./1e5, ':r');
    end

    % E-k for the valence bands ----------------------------------------
    figure(h_disp);
    k=k*10; %per Angstrom to per nm
    hold on;
    h=plot(k_t_vec,E1,'b');
    h=plot(k_t_vec,E2,'bx');
    h=plot(k_t_vec,E3,'b');
    h=plot(k_t_vec,E4,'b+');
    xlabel('k_t (m^-^1)');
    ylabel('Energy (eV)');
    grid on; box on;

    % Relative transition strengths ------------------------------------
    figure(h_rts);
    colors = [1 0 0; 1 0 1; 0 1 0; 0 1 1];
    for (jj=1:length(E_v{1}))
        l_tran1{jj} = sprintf('c1-h%d',jj);
        l_tran2{jj} = sprintf('c2-h%d',jj);

        subplot(221); hold on; grid on; box on;
        h_tran1(jj) = plot(k_t_vec, squeeze(rtsTE(:,1,jj)), 'Color', colors(jj,:));
        %h_tran1(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTE(:,1,jj)),30,'moving'), 'Color', colors(jj,:));
        %text(k_t_vec(end)*1e10, rtsTE(end,1,jj), sprintf('e1-h%d',jj));
        subplot(222); hold on; grid on; box on;
        h_tran2(jj) = plot(k_t_vec, squeeze(rtsTM(:,1,jj)), 'Color', colors(jj,:));
        %h_tran2(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTM(:,1,jj)),30,'moving'), 'Color', colors(jj,:));
        %text(k_t_vec(end)*1e10, rtsTM(end,1,jj), sprintf('e1-h%d',jj));
        subplot(223); hold on; grid on; box on;
        h_tran3(jj) = plot(k_t_vec, squeeze(rtsTE(:,2,jj)), 'Color', colors(jj,:));
        %h_tran3(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTE(:,2,jj)),30,'moving'), 'Color', colors(jj,:));
        %text(k_t_vec(end)*1e10, rtsTE(end,2,jj), sprintf('e2-h%d',jj));
        subplot(224); hold on; grid on; box on;
        h_tran4(jj) = plot(k_t_vec, squeeze(rtsTM(:,2,jj)), 'Color', colors(jj,:));
        %h_tran4(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTM(:,2,jj)),30,'moving'), 'Color', colors(jj,:));
        %plot(k_t_vec*1e10, spline(k_t_vec*1e10, squeeze(rtsTM(:,2,jj))));
        %text(k_t_vec(end)*1e10, rtsTM(end,2,jj), sprintf('e2-h%d',jj));
    end
    subplot(221); legend(h_tran1, l_tran1); title('TE'); xlabel('k_t (m^-^1)')
    subplot(222); legend(h_tran2, l_tran1); title('TM'); xlabel('k_t (m^-^1)')
    subplot(223); legend(h_tran3, l_tran2); title('TE'); xlabel('k_t (m^-^1)')
    subplot(224); legend(h_tran4, l_tran2); title('TM'); xlabel('k_t (m^-^1)');

    % Density of states -----------------------------------------------
    figure(h_dos); box on; grid on;
    Dos_t = zeros(1,length(Dos_E(1,:)));
    for (vb_index=1:length(E_v{1}))
        hold on;
        h_dos_fig(vb_index) = plot(E_grid*1e3, Dos_E(vb_index,:), 'Color', colors(vb_index,:));
        Dos_t = Dos_t + Dos_E(vb_index,:);
        hold off;
        l_dos{vb_index} = sprintf('h%d',vb_index);
    end
    hold on; 
    plot(E_grid*1e3, Dos_t,':','linewidth',3);
    hold off;
    xlabel('E (meV)'); ylabel('DOS (m^{-2}eV^{-1})');
    title('Valence band density of states (DOS)');
    legend(h_dos_fig,l_dos);

    % Gain ------------------------------------------------------------
    figure(h_gain); box on; grid on;
    ii = 1;
    for (vb_index=1:length(E_v{1}))
        for (cb_index=1:length(E_c))
            subplot(length(E_v{1}), length(E_c), ii);
            plot(E_pump/q, squeeze(Gain_E{vb_index, cb_index, :}));
            title(['cb', num2str(cb_index), '-vb', num2str(vb_index), ' - TE']);
            xlabel('E (meV)');  ylabel('Gain (cm^-^1)');
            ii = ii+1;
        end
    end

    % Absorption ------------------------------------------------------
    figure(h_abs); box on; grid on;
    hold on; 
    plot(1e3*E_grid/q, alpha_TE);
    plot(1e3*E_grid/q, alpha_TM, 'r');
    hold off;
    xlabel('E-E_g (meV)'); ylabel('\alpha');
    title('Absorption coefficient');
    
end