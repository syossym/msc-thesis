function Results = SingleQW_kp(well_width,x)

%
% This script calculates the band structure of a single QW structure, using
% the Luttinger-Kohn method.
%
%   Input: 'well_width' (A) - the width of the QW
%          'x' - Al concentration in the well AlGaAs alloy
%
%   Output: 'Results' - structure containing the major results of the
%                       simulation
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

run_cell = -1;    % marks the number of the cell to run. 
                  % Choose '-1' for the execution of the entire script.

                  
%% 1. Constants and definitions

if (run_cell == 1 || run_cell == -1)
   
   % Load the physical constants
    global Consts;

    % Creating the figures for plotting
    h_struct = figure('Name','Structure');
    h_env_func = figure('Name','Envelope Functions');
    h_disp = figure('Name', 'Dispersion');
    h_rts = figure('Name', 'Transition Matrix');

    % Simulation parameter definitions
    params.T = 2; params.E = 0;
    Mat_GaAs = GetMaterial('GaAs', params);

    Gw1 = (Consts.hbar^2)*Mat_GaAs.g1/(2*Consts.m_0);
    Gw2 = (Consts.hbar^2)*Mat_GaAs.g2/(2*Consts.m_0);
    Gw3 = (Consts.hbar^2)*Mat_GaAs.g3/(2*Consts.m_0);
    Cw = Gw1+Gw2;
    Dw = -2*Gw2+Gw1;
    Aw = Gw1-Gw2;
    Bw = 2*Gw2+Gw1;
    Ew = Gw2+Gw3;
    
    Mat_AlAs = GetMaterial('AlAs', params);

    a1 = (Consts.hbar^2)*Mat_AlAs.g1/(2*Consts.m_0);
    Gb1 = ((1-x)*Gw1)+(x*a1);
    a2 = (Consts.hbar^2)*Mat_AlAs.g2/(2*Consts.m_0);
    Gb2 = ((1-x)*Gw2)+(x*a2);
    a3 = (Consts.hbar^2)*Mat_AlAs.g3/(2*Consts.m_0);
    Gb3 = ((1-x)*Gw3)+(x*a3);
    Cb = Gb1+Gb2;
    Db = -2*Gb2+Gb1;
    Ab = Gb1-Gb2;
    Bb = 2*Gb2+Gb1;
    Eb = Gb2+Gb3;

    Ev=0; Evb=(((1-x)*0)+(x*0.75))*Consts.e_0;
    Eg = Mat_GaAs.E_g;   % eV

end

%% 2. Structure parameters

if (run_cell == 2 || run_cell == -1)

    dz = Mat_GaAs.a;         % here we consider the basic grid of the simulation
                             % to be the lattice constant
                             
    Nw = round(well_width/(dz*1e10) + 1);
                             
    %Nw = 30;
    Nb = Nw/2;
    Np = Nb+Nw+Nb;
    %W = (Nw-1)*dz*1e9;
    W = well_width;
    Z=zeros(Np);

    z = linspace(0,Np*dz,Np);
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

%% 3. Conduction band (using the shooting technique)

if (run_cell == 3 || run_cell == -1)

    hs_structure = [Nb*dz, Nb, x ; (W/1e10), Nw, 0 ; Nb*dz, Nb, x];
    [E_c,V_cb,V_vb,m_eff,con_profile] = ShootingMethod_E('e', hs_structure, 0, 2, Np, dz);
    wf_c  = ShootingMethod_WFs('e', E_c, V_cb, m_eff, con_profile);
    wf_c  = wf_c./sqrt(trapz(z, wf_c(1,:).^2));
    Evb = max(V_vb(:,2))./Consts.e_0;
    V_vb = -V_vb(:,2)./Consts.e_0;

    % Plot the structure -----------------------------------------------
    figure(h_struct);
    plot(z/1e-10, Eg*ones(1,Np) + V_cb(:,2)'./Consts.e_0, 'g', z/1e-10, V_vb, 'g', 'linewidth', 2);

end

%% 4. Velence band (4X4 k.p)

if (run_cell == 4 || run_cell == -1)

    for nk=1:100
        k(nk)=(nk-1)/500; % in A^-1
        
        l=0; m=1;
        lm=sqrt((l^2)+(m^2));
        kx=(l/lm)*k(nk)*1e10;
        ky=(m/lm)*k(nk)*1e10;
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
        En = diag(En)./Consts.e_0;
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

%% 5. Plotting  

if (run_cell == 5 || run_cell == -1)

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
        plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np)*Consts.e_0 + Fhh{1}(:,vb_index).'./1e5, ':b');
        plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np)*Consts.e_0 + Flh{1}(:,vb_index).'./1e5, ':r');
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
    
end

%% 6. Saving Results

if (run_cell == 6 || run_cell == -1)

    Results.m_eff = m_eff;
    Results.Np = Np;
    Results.Nw = Nw;
    Results.W = W;
    Results.k_t_vec = k_t_vec;
    Results.V_cb = V_cb;
    Results.V_vb = V_vb;
    Results.wf_c = wf_c;
    Results.E_c = E_c;
    Results.E_v = E_v;
    Results.rtsTE = rtsTE;
    Results.rtsTM = rtsTM;
    Results.rtsTE_k = rtsTE_k;
    Results.rtsTM_k = rtsTM_k;

end
