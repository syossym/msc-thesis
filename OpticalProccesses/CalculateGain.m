function OutParams = CalculateGain(vb_index, cb_index, calc_ver, InParams)

run_cell = -1;    % marks the number of the cell to run.
% Choose '-1' for the execution of the entire script.

%% 1. Constants and definitions

if (run_cell == 1 || run_cell == -1)
    
    global Consts;
    
    N = InParams.N;
    Materials = InParams.Materials;
    W_tot = 0;
    for (mm=1:length(Materials))
        W_tot = W_tot + Materials{mm}.Width*1e-10*100;  % [cm]
    end
    W = InParams.Materials{2}.Width*1e-10;    % [m]
    E_p = InParams.Materials{2}.E_p;          % [eV]
    z_grid = InParams.z_grid;
    k_t_vec = InParams.k_t_mat;
    m_eff_e = InParams.m_e_profile(length(InParams.m_e_profile)/2);
    V_vb = InParams.v_e_profile;
    v_h_profile = InParams.v_h_profile;
    E_v = InParams.E_k_v;
    E_c = InParams.E_0_c(:,2);
    rtsTE = InParams.rtsTE;
    rtsTM = InParams.rtsTM;
    Dos_h_k = InParams.Dos_h_k;
    Dos_h_E = InParams.Dos_h_E;
    
    E_v_k = zeros(length(E_v), length(k_t_vec));
    for (ii=1:length(E_v))
        for (nk = 1:length(k_t_vec))
            E_v_k(ii, nk) = abs(E_v{ii}(nk));   % J
        end
    end
    
    N_v = InParams.N_v;           % valence band population (1e10 cm^-2)
    N_c = InParams.N_c;           % conduction band population (1e10 cm^-2)
    Gamma = InParams.Gamma;       % Lorentzian linewidth
    T = InParams.T;               % temperature (K)
    
    meV = 1e-3*Consts.e_0;
    
    % Creating the plots
    h_calc = figure('Name', 'Calculation');
    h_gain = figure('Name', 'Gain');
    
end

%% 2. Calculation

if (run_cell == 2 || run_cell == -1)
    
    switch (calc_ver)
        
        case 1,
            
            E_grid_step = 0.01*1e-3*Consts.e_0;                    % [J]
            E_pump = [1.0*Consts.e_0:E_grid_step:2*Consts.e_0];    % grid for the pump energy
            % in excess of band gap energy (J)
            alpha_E_TE = zeros(length(E_pump), 1);
            alpha_E_TM = zeros(length(E_pump), 1);
            
            C_0 = (pi*Consts.e_0^2*Consts.hbar)/(3.3*Consts.c*Consts.eps_0*Consts.m_0^2);
            Mb_s = Consts.e_0*E_p*Consts.m_0/6;
            M_TE = Mb_s*squeeze(abs(rtsTE(vb_index,cb_index,:)));
            M_TM = Mb_s*squeeze(abs(rtsTM(vb_index,cb_index,:)));
            E_v_curr = -squeeze(E_v_k(vb_index,:));
            E_c_curr = 1.52*Consts.e_0 + E_c(cb_index)*Consts.e_0 + Consts.hbar^2*k_t_vec.^2/(2*0.0665*Consts.m_0);
            del_E_curr = E_c_curr - E_v_curr;
            
            dos_v = Dos_h_E(vb_index, Dos_h_E(vb_index,:) ~= 0); dos_v = dos_v(1)/Consts.e_0;
            Ef_v = CalculateQuasiFermiLevel('h', -E_v_curr(1) , N, T);
            %Ef_v = QuasiFermiLevels(max(abs(V_vb)), -E_v_curr(1), dos_v, N_v*1e10*1e4, T, 'Iterative');
            dos_c = 1e-4*0.0665*Consts.m_0/(pi*Consts.hbar^2);                                         % conduction band density of states [cm^-2]
            Ef_c = 1.52*Consts.e_0 + E_c(cb_index)*Consts.e_0 + Consts.k_B*T*log(exp(1e10*N_c/(dos_c*Consts.k_B*T))-1);  % conduction subband QF level, relative to the band edge
            f_v = FermiDirac(Ef_v, E_v_curr, T);
            f_c = FermiDirac(Ef_c, E_c_curr, T);
            del_f = f_c - f_v;
            
            for (aa = 1:length(E_pump))
                
                E_curr = E_pump(aa);
                B = (Gamma/(2*pi))./((del_E_curr-E_curr).^2+0.25*Gamma^2);
                alpha_E_TE(aa) = 4*C_0*(1/E_curr)*trapz(k_t_vec, (k_t_vec/(2*pi)).*M_TE'.*del_f.*B )/(pi*W);
                alpha_E_TM(aa) = 4*C_0*(1/E_curr)*trapz(k_t_vec, (k_t_vec/(2*pi)).*M_TM'.*del_f.*B )/(pi*W);
                
            end
            
        case 2,
            
            E_grid_step = 0.01*meV;                  % J
            E_pump = [0*meV:E_grid_step:200*meV];       % grid for the pump energy
            
            % in excess of band gap energy (J)
            gain_E_TE = zeros(length(E_pump), 1);
            gain_E_TM = zeros(length(E_pump), 1);
            
            [vb_index, cb_index]
            
            rts_TE = squeeze(rtsTE(vb_index,cb_index,:));
            rts_TM = squeeze(rtsTM(vb_index,cb_index,:));
            Dos_h_k_spec = squeeze(Dos_h_k(vb_index,:))/Consts.e_0;    % 1/(J m^2)
            E_v_curr = E_v_k(vb_index,:);
            E_c_curr = E_c(cb_index)*meV + Consts.hbar^2*k_t_vec.^2/(2*0.0665*Consts.m_0);
            
            gain_k_TE = (Consts.e_0^2)*Consts.hbar/(Consts.eps_0*Consts.c*Consts.m_0^2)*rts_TE'./(1./Dos_h_k_spec + pi*Consts.hbar^2/(0.0665*Consts.m_0));   % g(k)
            gain_k_TE = gain_k_TE./(1.52*Consts.e_0 + E_v_curr + E_c_curr);
            gain_k_TE = gain_k_TE*28.8*Consts.e_0*Consts.m_0/2;
            
            gain_k_TM = (Consts.e_0^2)*Consts.hbar/(Consts.eps_0*Consts.c*Consts.m_0^2)*rts_TM'./(1./Dos_h_k_spec + pi*Consts.hbar^2/(0.0665*Consts.m_0));   % g(k)
            gain_k_TM = gain_k_TM./(1.52*Consts.e_0 + E_v_curr + E_c_curr);
            gain_k_TM = gain_k_TM*28.8*Consts.e_0*Consts.m_0/2;

            % Calculating the quasi-Fermi energies
            DOS = Dos_h_E(vb_index, Dos_h_E(vb_index,:) ~= 0); DOS = DOS(1)/Consts.e_0;
            
            Ef_v = CalculateQuasiFermiLevel('h', E_v_curr(1), N_v, T, (DOS*pi*Consts.hbar^2));      % valence subband QF level, relative to the band edge
            Ef_c = E_c(cb_index)*meV + CalculateQuasiFermiLevel('e', E_c(cb_index)*meV, N_c, T);
            
            %Ef_v = QuasiFermiLevels(max(V_vb), E_n(1), DOS, N_v*1e10*1e4, T);
            %dos_c = 1e-4*0.0665*Consts.m_0/(pi*Consts.hbar^2);         % conduction band density of states [cm^-2]
            %Ef_c = E_c(cb_index)*meV + Consts.k_B*T*log(exp(1e10*N_c/(dos_c*Consts.k_B*T))-1);  % conduction subband QF level, relative to the band edge
            
            figure(h_calc);
            plot(k_t_vec, -E_v_curr/meV, 'b',...
                 k_t_vec, E_c_curr/meV, 'b',...
                 k_t_vec, -ones(size(k_t_vec)).*Ef_v/meV, ':g',...
                 k_t_vec, ones(size(k_t_vec)).*Ef_c/meV, ':r');
             
            f_coeff = (FermiDirac(Ef_v, E_v_curr, T) + ...
                       FermiDirac(Ef_c, E_c_curr, T) - ...
                       1);
            
            gain_k_TE = gain_k_TE.*f_coeff; gain_k_TM = gain_k_TM.*f_coeff;
            E_diff = E_v_k(vb_index,:) + E_c_curr;  % energy difference between the conduction and the valence bands vs. k
            dE_diff = diff(E_diff);
            
            % Fitting Gain(k) to Gain(E)
            for (step=1:size(dE_diff,2))
                E_1 = min(E_diff(step), E_diff(step+1));
                E_2 = max(E_diff(step), E_diff(step+1));
                
                node_index_1 = floor(E_1/E_grid_step+1);
                node_index_2 = ceil(E_2/E_grid_step+1);
                if (node_index_2 < length(E_pump))
                    for (node_index=node_index_1:(node_index_2-1))
                        gain_E_TE(node_index) = gain_E_TE(node_index) + gain_k_TE(step);
                        gain_E_TM(node_index) = gain_E_TM(node_index) + gain_k_TM(step);
                    end
                    node_fraction_1 = (E_1-E_pump(node_index_1))/E_grid_step;
                    gain_E_TE(node_index_1) = gain_E_TE(node_index_1) - node_fraction_1*gain_k_TE(step);
                    gain_E_TM(node_index_1) = gain_E_TM(node_index_1) - node_fraction_1*gain_k_TM(step);
                    node_fraction_2 = (E_pump(node_index_2)-E_2)/E_grid_step;
                    gain_E_TE(node_index_2-1) = gain_E_TE(node_index_2-1) - node_fraction_2*gain_k_TE(step); 
                    gain_E_TM(node_index_2-1) = gain_E_TM(node_index_2-1) - node_fraction_2*gain_k_TM(step);
                end
            end
            
            % Adding linewidth
            if (Gamma ~= 0)
                factor = 1/0.844;   % this factor makes up for the loss due to
                % limiting line shape width to 2*Gamma
                
                gain_E_lw_TE = zeros(length(E_pump), 1);
                gain_E_lw_TM = zeros(length(E_pump), 1);
                width = 2*Gamma;
                dE = E_grid_step;
                E_width = floor(width/dE);
                
                for (ii=1:length(E_pump))
                    for (jj=ii-E_width:ii+E_width)
                        if (jj>0 && jj<length(E_pump))
                            fract = dE*(Gamma/(2*pi))/((E_pump(ii)-E_pump(jj))^2+(Gamma/2)^2);
                            gain_E_lw_TE(jj) = gain_E_lw_TE(jj) + fract*gain_E_TE(ii);
                            gain_E_lw_TM(jj) = gain_E_lw_TM(jj) + fract*gain_E_TM(ii);
                        end
                    end
                end
                
                gain_E_TE = gain_E_lw_TE*factor;
                gain_E_TM = gain_E_lw_TM*factor;
            end
            
            alpha_E_TE = -gain_E_TE/W;
            alpha_E_TM = -gain_E_TM/W;
    end
end

%% 3. Plotting

if (run_cell == 3 || run_cell == -1)
    
    switch (calc_ver)
        
        case 1,
            
            figure(h_gain);
            plot(E_pump/(1e-3*Consts.e_0), 1e-2*alpha_E_TE, 'b', E_pump/(1e-3*Consts.e_0), 1e-2*alpha_E_TM, 'g');
            legend('TE', 'TM'); xlabel('E_{pump} [meV]'); ylabel('\alpha(E_{pump}) [cm^-^1]');
            title(['\alpha(E_{pump}) for [vb,cb] = [' num2str(vb_index) ',' num2str(cb_index) ']']);
            grid on; box on;
            
        case 2,
            
            figure(h_gain);
            plot(E_pump/(1e-3*Consts.e_0), 1e-2*alpha_E_TE, 'b', E_pump/(1e-3*Consts.e_0), 1e-2*alpha_E_TM, 'g');
            legend('TE', 'TM'); xlabel('E_{pump} [meV]'); ylabel('\alpha(E_{pump}) [cm^-^1]');
            title(['\alpha(E_{pump}) for [vb,cb] = [' num2str(vb_index) ',' num2str(cb_index) ']']);
            grid on; box on;
            
    end
end

%% 4. Saving Results

if (run_cell == 4 || run_cell == -1)
    
    OutParams.alpha_E_TE = alpha_E_TE;
    OutParams.alpha_E_TM = alpha_E_TM;
    OutParams.E_pump = E_pump;
    
end