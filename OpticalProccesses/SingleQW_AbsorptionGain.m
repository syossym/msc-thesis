function Results = SingleQW_AbsorptionGain(Struct)

%
% This script calculates the band structure of a single QW structure, using
% the Luttinger-Kohn method.
%
%   Input: 'Struct' - structure containing the QW parameters calculated by
%                     the SingleQW_kp function
%
%   Output: 'Results' - structure containing the major results of the
%                       simulation
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, October 2009
% Edited by: -
%

run_cell = -1;    % marks the number of the cell to run. 
                  % Choose '-1' for the execution of the entire script.
     
%% 1. Constants and definitions

if (run_cell == 1 || run_cell == -1)
   
   % Load the physical constants
    global Consts;
    
    Np = Struct.Np;
    m_eff = Struct.m_eff;
    k_t_vec = Struct.k_t_vec;
    V_cb = Struct.V_cb;
    V_vb = Struct.V_vb;
    wf_c = Struct.wf_c;
    E_c = Struct.E_c;
    E_v = Struct.E_v;
    rtsTE = Struct.rtsTE;
    rtsTM = Struct.rtsTM;
    rtsTE_k = Struct.rtsTE_k;
    rtsTM_k = Struct.rtsTM_k;
   
    Mat_GaAs = GetMaterial('GaAs');
    Eg = Mat_GaAs.E_g;
    
    % Creating the plots
    h_dos = figure('Name', 'DOS');
    h_abs = figure('Name', 'Absorption');
    h_gain = figure('Name', 'Gain');
    
end
                  
%% 2. DOS

if (run_cell == 2 || run_cell == -1)

    E_grid_step = 0.1e-3;                                   % eV
    k_t_step = k_t_vec(2)-k_t_vec(1);

    E_v_k = zeros(length(E_v{1}), length(k_t_vec));
    for (vb_index=1:length(E_v{1}))
        for (nk = 1:length(k_t_vec))
            E_v_k(vb_index, nk) = abs(E_v{nk}(vb_index));   % eV
        end
    end

    E_min = min(abs(V_vb));                          % eV
    E_max = max(abs(V_vb));                          % eV
    E_grid_DOS = [E_min:E_grid_step:E_max];              % energy grid for the DOS (eV)
    E_grid_size = size(E_grid_DOS,2);
    Dos_E = zeros(length(E_v{1}), length(E_grid_DOS));
    dE_k = diff(E_v_k,1,2);                          % eV
    dk_t = diff(k_t_vec);                            % m
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
                node_fraction_1 = (E_1-E_grid_DOS(node_index_1))/E_grid_step;
                Dos_E(vb_index,node_index_1) = Dos_E(vb_index,node_index_1) - node_fraction_1*Dos_k(vb_index,step);
                node_fraction_2 = (E_grid_DOS(node_index_2)-E_2)/E_grid_step;
                Dos_E(vb_index,node_index_2-1) = Dos_E(vb_index,node_index_2-1) - node_fraction_2*Dos_k(vb_index,step);
            end
        end
    end

    diff_Dos_k = diff(Dos_k,1,2);
    Dos_k = [Dos_k(:,1), Dos_k(:,1:length(k_t_vec)-2)+diff_Dos_k, Dos_k(:,length(k_t_vec)-1)];

end

%% 3. Gain 

if (run_cell == 3 || run_cell == -1)

    % Parameters -------------------------------------------------------
    N_v = 4e2;                              % valence band population (1e10 cm^-2)
    N_c = 5e2;                              % conduction band population (1e10 cm^-2)
    Gamma = 2e-3*Consts.e_0;                % Lorentzian linewidth
    T = 300;                                % temperature (K)

    E_grid_step = 0.1*1e-3*Consts.e_0;      % J

    E_pump = [0:E_grid_step:500*1e-3*Consts.e_0];  % grid for the pump energy 
                                                   % in excess of band gap energy (J)
    gain_E = zeros(length(E_pump), 1);

    % Gain calculation for all valence & conduction combinations -------
    for (vb_index = 1:length(E_v{1}))
        for (cb_index=1:length(E_c))
            
            [vb_index, cb_index]
            
            rts = squeeze(rtsTE_k(vb_index,cb_index,:));
            Dos_k_spec = squeeze(Dos_k(vb_index,:))/Consts.e_0;    % 1/(J m^2) 
            E_n = squeeze(E_v_k(vb_index,:)).*Consts.e_0;

            gain_k = (Consts.e_0^2)*Consts.hbar/(Consts.eps_0*Consts.c*Consts.m_0^2)*rts'./(1./Dos_k_spec + pi*Consts.hbar^2/(0.0665*Consts.m_0));   % g(k)
            gain_k = gain_k./(1.52*Consts.e_0 + E_n + E_c(cb_index)*Consts.e_0 + Consts.hbar^2*k_t_vec.^2/(2*0.0665*Consts.m_0));
            gain_k = gain_k*28.8*Consts.e_0*Consts.m_0/2;

            I = find(Dos_E(vb_index, :) ~= 0);
            dE_f_v = QuasiFermiLevels((V_vb)*Consts.e_0, (E_v{1}(vb_index))*Consts.e_0, Dos_E(vb_index, I(1))/Consts.e_0, 1e4*N_v*1e10, T);
            E_f_v = abs(E_v{1}(vb_index))*Consts.e_0 + dE_f_v;
            
            % QuasiFermiLevels(vb_index, N_v, T);         % valence subband quasi-Fermi level, relative to the valence band edge

            Dos_c = 1e-4*m_eff(Np/2,2)/(pi*Consts.e_0*Consts.hbar^2);                       % density of states in the conduction band (cm^-2 J^-1)
            dE_f_c = Consts.k_B_J*T*log(exp(1e10*N_c/(Dos_c*Consts.k_B_J*T))-1);    % conduction subband quasi-Fermi level, relative to the conduction subband edge
            E_f_c = E_c(cb_index)*Consts.e_0 + dE_f_c;                              % conduction band quasi-Fermi level, relative to the conduction band edge

            gain_k = gain_k.*(FermiDirac(E_f_v, E_v_k(vb_index,:)*Consts.e_0, T) + ...
                              FermiDirac(E_f_c, E_c(cb_index)*Consts.e_0+Consts.hbar^2*k_t_vec.^2/(2*0.0665*Consts.m_0), T) - ...
                              1);
            E_diff = E_v_k(vb_index,:)*Consts.e_0 + E_c(cb_index)*Consts.e_0 + Consts.hbar^2*k_t_vec.^2/(2*0.0665*Consts.e_0);  % energy difference between the conduction and the valence bands vs. k
            dE_diff = diff(E_diff);

            % Fitting Gain(k) to Gain(E)
            for (step=1:size(dE_diff,2))
                E_1 = min(E_diff(step), E_diff(step+1));
                E_2 = max(E_diff(step), E_diff(step+1));

                node_index_1 = floor(E_1/E_grid_step+1);
                node_index_2 = ceil(E_2/E_grid_step+1);
                if (node_index_2 < length(E_pump))
                    for (node_index=node_index_1:(node_index_2-1))
                        gain_E(node_index) = gain_E(node_index) + gain_k(step);
                    end
                    node_fraction_1 = (E_1-E_pump(node_index_1))/E_grid_step;
                    gain_E(node_index_1) = gain_E(node_index_1) - node_fraction_1*gain_k(step);
                    node_fraction_2 = (E_pump(node_index_2)-E_2)/E_grid_step;
                    gain_E(node_index_2-1) = gain_E(node_index_2-1) - node_fraction_2*gain_k(step);
                end
            end

            % Adding linewidth
            if (Gamma ~= 0)
                factor = 1/0.844;   % this factor makes up for the loss due to
                                    % limiting line shape width to 2*Gamma

                gain_E_lw = zeros(length(E_pump), 1);
                width = 2*Gamma;
                dE = E_grid_step;
                E_width = floor(width/dE);

                for (ii=1:length(E_pump))
                    for (jj=ii-E_width:ii+E_width)
                        if (jj>0 && jj<length(E_pump))
                            fract = dE*(Gamma/(2*pi))/((E_pump(ii)-E_pump(jj))^2+(Gamma/2)^2);
                            gain_E_lw(jj) = gain_E_lw(jj) + fract + gain_E(ii);
                        end
                    end
                end

                gain_E = gain_E_lw*factor;
            end

            Gain_E{vb_index, cb_index, :} = gain_E;
            Gain_k{vb_index, cb_index, :} = gain_k;
        end
    end
end

%% 4. Absorption 

if (run_cell == 4 || run_cell == -1)

    sums_TE = [];
    sums_TM = [];
    
    Gamma = 2e-3*Consts.e_0;       % J
    w_grid = k_t_vec*Consts.c;
    w_grid(1) = w_grid(2);
    E_grid_Abs = w_grid*Consts.hbar;   % J
    
    rts_TE = abs(squeeze(rtsTE_k(:,1,:)));
    rts_TM = abs(squeeze(rtsTM_k(:,1,:)));
    
    for (vb_index = 1:length(E_v{1}))   
        del_E = (E_c(1) + Eg + E_v{1}(vb_index))*Consts.e_0;    % J
        line_width = (Gamma/(2*pi))./((del_E*ones(size(E_grid_Abs))-E_grid_Abs).^2 ... 
                     + ones(size(E_grid_Abs))*(Gamma/2)^2);
        sums_TE = [sums_TE; rts_TE(vb_index,:).*line_width];
        sums_TM = [sums_TM; rts_TM(vb_index,:).*line_width];
    end

    alpha_TE = (1./w_grid) .* sum(sums_TE);
    alpha_TM = (1./w_grid) .* sum(sums_TM);

end


%% 5. Plotting

if (run_cell == 5 || run_cell == -1)

    colors = [1 0 0; 1 0 1; 0 1 0; 0 1 1];
    
    % Density of states -----------------------------------------------
    figure(h_dos); box on; grid on;
    Dos_t = zeros(1,length(Dos_E(1,:)));
    for (vb_index=1:length(E_v{1}))
        hold on;
        h_dos_fig(vb_index) = plot(E_grid_DOS*1e3, Dos_E(vb_index,:), 'Color', colors(vb_index,:));
        Dos_t = Dos_t + Dos_E(vb_index,:);
        hold off;
        l_dos{vb_index} = sprintf('h%d',vb_index);
    end
    hold on;
    plot(E_grid_DOS*1e3, Dos_t,':','linewidth',3);
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
            plot(E_pump/Consts.e_0, squeeze(Gain_E{vb_index, cb_index, :}));
            title(['cb', num2str(cb_index), '-vb', num2str(vb_index), ' - TE']);
            xlabel('E (meV)');  ylabel('Gain (cm^-^1)');
            ii = ii+1;
        end
    end

    % Absorption ------------------------------------------------------
    figure(h_abs); box on; grid on;
    hold on;
    plot(1e3*E_grid_Abs/Consts.e_0, alpha_TE);
    plot(1e3*E_grid_Abs/Consts.e_0, alpha_TM, 'r');
    hold off;
    xlabel('E-E_g (meV)'); ylabel('\alpha');
    title('Absorption coefficient');
    
    
end

%% 6. Saving Results

if (run_cell == 6 || run_cell == -1)

    Results.Dos_k = Dos_k;
    Results.Dos_E = Dos_E;
    Results.Gamma = Gamma;
    Results.alpha_TE = alpha_TE;
    Results.alpha_TM = alpha_TM;
    Results.Gain_E = Gain_E;
    Results.Gain_k = Gain_k;
    Results.E_grid_DOS = E_grid_DOS;
    Results.E_grid_Abs = E_grid_Abs;
    Results.E_pump = E_pump; 
    
end


