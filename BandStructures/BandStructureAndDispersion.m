function Results = BandStructureAndDispersion(Structure, method)

%
% This script calculates the band structure of an arbitrary Quantum
% structure using a 4X4 (diagonalized) k.p and the propagation matrix method. The shooting
% method and other codes written in C are used throughout the simulation.
%
%   Input: 'Structure' - cell array containing the detailes of the
%                        simulated structure.
%
%   Output: 'Results' - struct containing the major results of the
%                       simulation.
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
    
    % Init the material parameters structure
    [row_num,column_num] = size(Structure);
    Materials = cell(row_num,1);
    if (column_num == 2)
        for (ii=1:row_num)
            Materials{ii} = GetMaterial(Structure{ii,1});
            Materials{ii}.Width = Structure{ii,2};
        end
    elseif (column_num == 3)
        for (ii=1:length(Structure(:,1)))
            mat_params.x = Structure{ii,3};
            Materials{ii} = GetMaterial(Structure{ii,1},mat_params);
            Materials{ii}.Width = Structure{ii,2};
        end
    else
        error('TransferMatrix:Input_parameter_error',...
            'The input cell array Structure must have 2 or 3 columns');
    end
    
    % Create a 's.r' file for the structure
    delete('.\Bin\Input\s.r');
    delete('.\Bin\*.r');
    s_file = zeros(row_num,3);
    for (ii=1:row_num)
        s_file(ii,:) = [Materials{ii}.Width, Materials{ii}.x, 0];
    end
    save '.\Bin\Input\s.r' 's_file' -ASCII;
    copyfile('.\Bin\Input\s.r', '.\Bin\s.r');
    
    % Generate the structure profiles file 'x.r'
    RunExternalCode('.\Bin\efsx.exe');
    copyfile('.\Bin\x.r', '.\Bin\Output\x.r');
    
    % Generate the potential profiles 'v.r'
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxv.exe', ' -p e');
    copyfile('.\Bin\v.r', '.\Bin\Output\v_e.r');
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxv.exe', ' -p h');
    copyfile('.\Bin\v.r', '.\Bin\Output\v_h.r');
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxv.exe', ' -p l');
    copyfile('.\Bin\v.r', '.\Bin\Output\v_l.r');
    
    % Generate effective mass profiles 'm.r'
    delete('.\Bin\m.r');
    RunExternalCode('.\Bin\efxm.exe', ' -p e');
    copyfile('.\Bin\m.r', '.\Bin\Output\m_e.r');
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxm.exe', ' -p h');
    copyfile('.\Bin\m.r', '.\Bin\Output\m_h.r');
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxm.exe', ' -p l');
    copyfile('.\Bin\m.r', '.\Bin\Output\m_l.r');
    
    % Calculate the electronic levels and WF's using the shooting method
    delete('.\Bin\v.r');
    RunExternalCode('.\Bin\efxv.exe', ' -p e');
    delete('.\Bin\m.r');
    RunExternalCode('.\Bin\efxm.exe', ' -p e');
    delete('.\Bin\Ee.r');
    RunExternalCode(['.\Bin\efshoot.exe', ' -p e -s ',...
                     num2str(6)]);    % calculating levels
    copyfile('.\Bin\Ee.r', '.\Bin\Output\Ee.r');
    RunExternalCode('.\Bin\efwf.exe');             % calculating wf's
    copyfile('.\Bin\wf_e*', '.\Bin\Output\');
    
    % Init internal data structures
    xr_file = load('.\Bin\Output\x.r');
    z_grid = xr_file(:,1);                         % [m]
    x_profile = xr_file(:,2);
    vr_e_file = load('.\Bin\Output\v_e.r');
    v_e_profile = vr_e_file(:,2);                  % [eV]
    vr_h_file = load('.\Bin\Output\v_h.r');
    v_h_profile = vr_h_file(:,2);                  % [eV]
    mr_e_file = load('.\Bin\Output\m_e.r');
    m_e_profile = mr_e_file(:,2);                  % [kg]
    mr_h_file = load('.\Bin\Output\m_h.r');
    m_h_profile = mr_h_file(:,2);                  % [kg]
    E_0_c = load('.\Bin\Output\Ee.r');             % [meV]
    for (cb_num = 1:6)
        eval(['wf_e_temp=load(''.\Bin\Output\wf_e' num2str(cb_num) '.r'');']);
        wf_e_mat{cb_num} = wf_e_temp;
    end
    
    g1_profile = zeros(length(z_grid),1);
    g2_profile = zeros(length(z_grid),1);
    g3_profile = zeros(length(z_grid),1);
    g1_profile(1) = Materials{1}.g1;
    g2_profile(1) = Materials{1}.g2;
    g3_profile(1) = Materials{1}.g3;
    
    ii = 2; mat_index = 1;
    while(ii <= length(z_grid) && mat_index <= row_num)
        g1_profile(ii) = Materials{mat_index}.g1;
        g2_profile(ii) = Materials{mat_index}.g2;
        g3_profile(ii) = Materials{mat_index}.g3;
        if (x_profile(ii)~=x_profile(ii-1))
            mat_index=mat_index+1;
        end
        ii=ii+1;
    end
    
    norm_const = Consts.hbar^2/(2*Consts.m_0);
    gamma1_profile = g1_profile.*norm_const;
    gamma2_profile = g2_profile.*norm_const;
    gamma3_profile = g3_profile.*norm_const;
    
    Profile.z_grid = z_grid;
    Profile.x_profile = x_profile;
    Profile.v_e_profile = v_e_profile;
    Profile.v_h_profile = v_h_profile;
    Profile.m_e_profile = m_e_profile;
    Profile.m_h_profile = m_h_profile;
    Profile.g1_profile = g1_profile;
    Profile.g2_profile = g2_profile;
    Profile.g3_profile = g3_profile;
    Profile.gamma1_profile = gamma1_profile;
    Profile.gamma2_profile = gamma2_profile;
    Profile.gamma3_profile = gamma3_profile;
    
end

%% 2. Simulation parameters

if (run_cell == 2 || run_cell == -1)
    
    % Common parameters
    meV = 1e-3*Consts.e_0;         % [J]
    dz = z_grid(2)-z_grid(1);      % grid step size [m]
    dzz = dz^2;
    N  = size(v_h_profile,1);      % num of grid points  
    k_t_max = 1e9;                 % max. calculated k_t
    k_t_step = 2e6;                % k_t step for E-k dispersion
    k_t_mat = [(0:0.25:1.75),...
               2:(k_t_max/k_t_step)]*k_t_step;
    num_cb_subbands = 2;           % number of conduction subbands                                  
    num_vb_subbands = 4;           % number of valence subbands    
    
    % Method specific parameters
    switch(method)
        case 'Matrix'
            
            %dz = z_grid(2)-z_grid(1);      % grid step size [m]
            %N  = size(v_h_profile,1);      % num of grid points
            nrm = 1e7;                     % normalization coefficient
                        
            min_E_diff = 0.05*meV;         % helps to determine whether two solutions are the same [J]
            a_lim = 1e5;                   % max. for a
            k_t = 1e5;                     % k_t for band-edge calculation
            E_min = min(v_h_profile);      % energy minimum [J]
            E_max = max(v_h_profile); % min(v_h_profile(1),...
            %v_h_profile(N));   % energy maximum [J]
            E_step = 0.1*meV;    	       % energy step for band-edge calculation [J]
            E_acc  = 1e-9*meV;             % energy resolution of the calculated subband energies [J]
            
            % Search parameters
            subband_E_search  = 0.15*meV;  % energy search range [J]
            E_abort           = meV;       % abort search when E range > E_abort
            %k_t_max           = 1e9;       % max. calculated k_t
            %k_t_step          = 2e6;       % k_t step for E-k dispersion
            %num_cb_subbands   = 2;         % number of conduction band subbands taken into account
            %num_vb_subbands   = nr-1;      % number of valence band subbands taken into account
            
            % Plotting parameters
            cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';...
                'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];
            ch = cl;
            
        case '4X4_k.p'
            
            Np = length(v_h_profile);
            A = gamma1_profile - gamma2_profile;
            C = gamma1_profile + gamma2_profile;
            B = 2*gamma2_profile + gamma1_profile;
            D = -2*gamma2_profile + gamma1_profile;
            E = gamma2_profile + gamma3_profile;
            G2 = gamma2_profile;
            G3 = gamma3_profile;
            V = v_h_profile;
            
        otherwise
            disp('Unrecognized calculation method');
    end
    
end

%% 3. Plotting definitions

if (run_cell == 3 || run_cell == -1)
    
    switch(method)
        case 'Matrix'
            
            h_LH_Subbnads = figure('Name','LH Subbnads');
            h_HH_Subbnads = figure('Name','HH Subbnads');
            h_E_k = figure('Name','E_k');
            h_rts_TE = figure('Name','RTS TE');
            h_rts_TM = figure('Name','RTS TM');
            
        case '4X4_k.p'
            
            h_struct = figure('Name','Structure');
            h_env_func = figure('Name','Envelope Functions');
            h_disp = figure('Name', 'Dispersion');
            h_rts = figure('Name', 'Transition Matrix');
            
        otherwise
            disp('Unrecognized calculation method');
    end
    
end

%% 4. Valence band-edge subband energy & WF calculation

if (run_cell == 4 || run_cell == -1)
    
    switch(method)
        case 'Matrix'
            
            % Prepare the plots
            figure(h_LH_Subbnads);
            plot(z_grid/1e-10, v_h_profile/meV, 'k');
            hold; title('LH Subbands'); xlabel('z [A]'); ylabel('E [meV]');
            figure(h_HH_Subbnads);
            plot(z_grid/1e-10, v_h_profile/meV, 'k');
            hold; title('HH Subbands'); xlabel('z [A]'); ylabel('E [meV]');
            
            % Start energy loop
            nr = 1;
            old_sign_l = 1;
            old_sign_h = 1;
            old_sign_d = 1;
            old_sign_d2 = 1;
            old_sign_d3 = 1;
            ii = 1;
            
            for (E = E_min:E_step:E_max)
                %disp(['Step ' num2str(ii) ': E=' num2str(E)]);
                
                % Transfer matrix
                tf_mat = TransferMatrix(E, k_t, Profile);
                h_10(ii) = tf_mat(3,3);
                l_10(ii) = tf_mat(4,3);
                h_01(ii) = tf_mat(3,4);
                l_01(ii) = tf_mat(4,4);
                c_min(ii) = -(tf_mat(1,3)*tf_mat(1,4) + tf_mat(2,3)*tf_mat(2,4)) / (tf_mat(1,4)^2 + tf_mat(2,4)^2);
                diff_h_min(ii) = tf_mat(3,3) + c_min(ii)*tf_mat(3,4);
                diff_l_min(ii) = tf_mat(4,3) + c_min(ii)*tf_mat(4,4);
                d_min(ii) = diff_h_min(ii)^2 + diff_l_min(ii)^2;
                dum(ii) = abs(h_10(ii)*h_01(ii) + l_01(ii)*l_10(ii));
                
                if (ii>1)
                    diff_d_min(ii)  = d_min(ii) - d_min(ii-1);
                    new_sign_d      = sign(diff_d_min(ii));
                    
                    if (ii>2)
                        if (new_sign_d-old_sign_d==2)
                            
                            [a, E_result] = WavefunctionMinimum(Profile,E-2*E_step,E,k_t,E_acc);
                            
                            %disp(['E_result=' num2str(E_result/meV)]);
                            
                            if (nr==1) | ( (nr>1) & abs(E_result-E_state(nr-1))>min_E_diff )
                                E_state(nr) = E_result;
                                [f_h, f_l]  = PlotWavefunction(E_result,k_t,1,a,Profile);
                                
                                %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);
                                
                                if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4 )
                                    E_state(nr) = E_result;
                                    coeff(nr)   = a;
                                    wf_l(nr,:)  = f_l;
                                    wf_h(nr,:)  = f_h;
                                    
                                    % Plotting
                                    figure(h_LH_Subbnads);
                                    plot(z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                                    drawnow;
                                    refresh(h_LH_Subbnads);
                                    figure(h_HH_Subbnads);
                                    plot(z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                                    drawnow;
                                    refresh(h_HH_Subbnads);
                                    
                                    wf_h_mat{nr} = [z_grid, f_h', f_l'];
                                    
                                    nr = nr + 1;
                                end
                            end
                        end
                    end
                    
                    old_sign_d = new_sign_d;
                end
                
                if (ii>1)
                    diff3(ii)   = dum(ii) - dum(ii-1);
                    new_sign_d3 = sign(diff3(ii));
                    if (ii>2)
                        if (new_sign_d3-old_sign_d3==2)
                            
                            [a, E_result] = WavefunctionMinimum(Profile,E-2*E_step,E+2*E_step,k_t,E_acc,E_step/50);
                            
                            if (nr==1) | ( (nr>1) & abs(E_result-E_state(nr-1))>min_E_diff )
                                
                                E_state(nr) = E_result;
                                [f_h,f_l]   = PlotWavefunction(E_result,k_t,1,a,Profile);
                                
                                %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);
                                
                                if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4  )
                                    E_state(nr) = E_result;
                                    coeff(nr)   = a;
                                    wf_l(nr,:)  = Normalize(f_l);
                                    wf_h(nr,:)  = Normalize(f_h);
                                    
                                    % Plotting
                                    figure(h_LH_Subbnads);
                                    plot(z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                                    drawnow;
                                    refresh(h_LH_Subbnads);
                                    figure(h_HH_Subbnads);
                                    plot(z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                                    drawnow;
                                    refresh(h_HH_Subbnads);
                                    
                                    wf_h_mat{nr} = [z_grid, f_h', f_l'];
                                    
                                    nr = nr + 1;
                                    
                                else
                                    
                                    [a, E_result] = WavefunctionMinimum(Profile,E_result,E+2*E_step,k_t,E_acc,E_step/50,2);
                                    E_state(nr) = E_result;
                                    [f_h,f_l]   =  PlotWavefunction(E_result,k_t,1,a,Profile);
                                    
                                    %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);
                                    
                                    if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4 )
                                        E_state(nr) = E_result;
                                        wf_l(nr,:)  = Normalize(f_l);
                                        wf_h(nr,:)  = Normalize(f_h);
                                        
                                        % Plotting
                                        figure(h_LH_Subbnads);
                                        plot(z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                                        drawnow;
                                        refresh(h_LH_Subbnads);
                                        figure(h_HH_Subbnads);
                                        plot(z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                                        drawnow;
                                        refresh(h_HH_Subbnads);
                                        
                                        wf_h_mat{nr} = [z_grid, f_h', f_l'];
                                        
                                        nr = nr + 1;
                                    end
                                end
                            end
                        end
                    end
                    
                    old_sign_d3 = new_sign_d3;
                end
                
                E_index(ii) = E;
                ii = ii+1;
            end
            
        case '4X4_k.p'
            
            
            
        otherwise
            disp('Unrecognized calculation method');
    end
    
end

%% 5. E-k dispersion relation calculation

if (run_cell == 5 || run_cell == -1)
    
    switch(method)
        case 'Matrix'
            
            num_vb_subbands = nr - 1;
            k_t_len = [];
            
            % Valence band calculations
            for (nr = 1:num_vb_subbands)
                clear  subband_E  coeff;
                
                % Definitions
                E_0 = E_state(nr);                    % subband edge energy
                subband_E(1) = E_0;
                sbnd_E_srch_act = subband_E_search;   % search range
                sbnd_E_srch_old = sbnd_E_srch_act;    % old search range (previous k)
                E_step = subband_E_search/10;         % energy step in range
                E_step_old = E_step;                  % old energy step
                dE = 0;                               % difference between two previous energies
                jj = 1;
                abort = 0;
                err = 0;                              % 1 if no minimum found
                sbnd_E_lim = meV*3e-6;                % minimum search range
                old_err = 0;                          % diff between the predicted and calculated E
                accuracy_modifier = 1;
                min_num = 1;                          % number of minima in search range
                E_acc_act = E_acc;
                %k_t_mat = [(0:0.25:1.75)  2:(k_t_max/k_t_step)]*k_t_step;
                
                tic;
                
                % Iterating over k_t
                while ((jj < size(k_t_mat,2)) && ~abort)
                    
                    k_t = k_t_mat(jj+1);
                    
                    min_E = subband_E(jj) + dE - sbnd_E_srch_act;
                    max_E = subband_E(jj) + dE + sbnd_E_srch_act;
                    min_E_arc(jj+1) = min_E;
                    max_E_arc(jj+1) = max_E;
                    
                    [a,E_result,err] = WavefunctionMinimum(Profile,min_E,max_E,k_t,E_acc_act,E_step,min_num);
                    
                    if (jj/50==floor(jj/50))
                        drawnow;
                    end
                    
                    if (err)    % -- minimum not found
                        
                        if ((min_num==2) && (jj<4))
                            min_num = 1;
                        else
                            num_steps = 2*sbnd_E_srch_act/E_step;
                            
                            if ((E_step>E_step_old/10) && (sbnd_E_srch_act<sbnd_E_srch_old*10))
                                
                                if (num_steps<3e2)
                                    E_step = E_step/2;
                                    num_steps = sbnd_E_srch_act/E_step;
                                    %disp(['Adjusting E_step to ', num2str(E_step/meV), ', ', num2str(floor(num_steps)), ' steps']);
                                else
                                    sbnd_E_srch_act = sbnd_E_srch_act*2;
                                    E_step = sbnd_E_srch_act/10;
                                    %disp(['Too many steps - increaing search area to ', num2str(sbnd_E_srch_act/meV), ', E_step=', num2str(E_step/meV), ', 10 steps']);
                                end
                                
                            elseif (sbnd_E_srch_act<5*sbnd_E_srch_old)
                                
                                sbnd_E_srch_act = sbnd_E_srch_old*5;
                                E_step = E_step_old/2;
                                num_steps = sbnd_E_srch_act/E_step;
                                %disp(['Search area too small - increasing search area to ', num2str(sbnd_E_srch_act/meV), ', E_step=', num2str(E_step/meV), ', ', num2str(num_steps), ' steps']);
                                
                            else
                                
                                k_t_mat_size = size(k_t_mat, 2);
                                k_t_step_local = k_t - k_t_mat(jj);
                                k_t_inter = k_t - k_t_step_local/2;
                                k_t_mat = [k_t_mat(1:jj), k_t_inter, k_t_mat(jj+1:k_t_mat_size)];
                                
                                P = polyfit(k_t_mat((jj-3):jj)/k_t_step, subband_E((jj-3):jj)/meV, 3);
                                E_predict = polyval(P, k_t_mat(jj+1)/k_t_step)*meV;
                                dE = E_predict - subband_E(jj);
                                sbnd_E_srch_act = sbnd_E_srch_old;
                                E_step = E_step_old;
                                
                                %disp(['Inserting additional k_t=', num2str(k_t)]);
                            end
                            
                            if (sbnd_E_srch_act>0.25*meV)
                                k_t_mat_size = size(k_t_mat, 2);
                                k_t_step_local = k_t - k_t_mat(jj);
                                k_t_inter = k_t - k_t_step_local/2;
                                k_t_mat = [k_t_mat(1:jj), k_t_inter, k_t_mat(jj+1:k_t_mat_size)];
                                P = polyfit(k_t_mat((jj-3):jj)/k_t_step, subband_E((jj-3):jj)/meV, 3);
                                E_predict = polyval(P, k_t_mat(jj+1)/k_t_step)*meV;
                                dE = E_predict - subband_E(jj);
                                sbnd_E_srch_act = sbnd_E_srch_old;
                                E_step = E_step_old;
                                %disp(['Search area too large, inserting additional k_t=', num2str(k_t)]);
                            end
                            
                            err = 0;
                        end
                    else                        % -- minimum found
                        
                        [f_h,f_l] = PlotWavefunction(E_result, k_t, 1, a, Profile);
                        
                        if ( ((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2)) < 1e-4*accuracy_modifier )
                            subband_E(jj+1) = E_result;
                            coeff(jj+1) = a;
                            
                            if (jj<4)
                                dE = (E_result - subband_E(jj))*(k_t_mat(jj+2)-k_t)/(k_t-k_t_mat(jj));
                            elseif (jj<(size(k_t_mat,2)-1))
                                P = polyfit(k_t_mat((jj-2):(jj+1))/k_t_step, subband_E((jj-2):(jj+1))/meV, 3);
                                E_predict = polyval(P, k_t_mat(jj+2)/k_t_step)*meV;
                                dE = E_predict - E_result;
                            end
                            
                            for (cb_num=1:num_cb_subbands)
                                overlap1 = trapz(z_grid, (wf_e_mat{cb_num}(:,2).*f_h'));
                                overlap2 = trapz(z_grid, (wf_e_mat{cb_num}(:,2).*f_l'));
                                
                                rtsTE(nr, cb_num, jj+1) = (3/4)*((abs(overlap1^2) + (1/3)*abs(overlap2^2)));
                                rtsTM(nr, cb_num, jj+1) = abs(overlap2^2);
                            end
                            
                            err = E_result - min_E - sbnd_E_srch_act;
                            jj = jj + 1;
                            
                            %                     disp(['k_t=', num2str(k_t), ', E_acc=', num2str(E_acc_act/meV), ...
                            %                           ', srch_dE=', num2str(sbnd_E_srch_act/meV), ' E_step=', num2str(E_step/meV) ...
                            %                           ' E found ', num2str(E_result/meV), ' meV']);
                            
                            sbnd_E_srch_old = sbnd_E_srch_act;
                            sbnd_E_srch_act = min([max(sbnd_E_srch_act/5, abs(5*err)), subband_E_search]);
                            sbnd_E_srch_fnd = sbnd_E_srch_act;
                            
                            if (sbnd_E_srch_act < sbnd_E_lim)
                                sbnd_E_srch_act = sbnd_E_lim;
                            end
                            
                            sbnd_E_lim = max([sbnd_E_lim/10, abs(5*err), min(E_acc,E_acc_act*1e4)]);
                            
                            E_step_old = E_step;
                            E_step = sbnd_E_srch_act/10;
                            E_step_fnd = E_step;
                            
                            E_acc_prev = E_acc_act;
                            E_acc_act = min(E_acc*1e5, E_acc_act*1.25);
                            E_acc_fnd = E_acc_act;
                            
                            accuracy_modifier = max(accuracy_modifier*0.75, 1);
                            changed_num = 0;
                            
                        else
                            
                            % Adjusting search parameters
                            if (jj<3)
                                E_acc_act = E_acc_act/10;
                                
                                if (E_acc_act<E_acc/20)
                                    if (min_num==1)
                                        min_num = 2;
                                    else
                                        min_num = 1;
                                    end
                                end
                            else
                                E_acc_act = E_acc_act*0.5;
                                %disp(['E_acc changed to ', num2str(E_acc_act/meV)]);
                                accuracy_modifier = min(accuracy_modifier*2, 100);
                                %disp(['accuracy_modifier=', num2str(accuracy_modifier)]);
                                
                                if (E_acc_act<E_acc_prev/3)
                                    
                                    k_t_mat_size = size(k_t_mat, 2);
                                    k_t_step_local = k_t - k_t_mat(jj);
                                    k_t_inter = k_t - k_t_step_local/2;
                                    k_t_mat = [k_t_mat(1:jj), k_t_inter, k_t_mat(jj+1:k_t_mat_size)];
                                    
                                    if (jj>=3)
                                        if (jj==3)
                                            index = 1;
                                        else
                                            index = jj-3;
                                        end
                                        P = polyfit(k_t_mat(index:jj)/k_t_step, subband_E(index:jj)/meV, 3);
                                        E_predict = polyval(P, k_t_mat(jj+1)/k_t_step)*meV;
                                    end
                                    
                                    dE = E_predict - subband_E(jj);
                                    sbnd_E_srch_act = sbnd_E_srch_old;
                                    E_step = E_step_old;
                                    
                                    %disp(['Inserting additional k_t=', num2str(k_t_inter)]);
                                end
                                
                                if (E_acc_act<E_result*eps*10)
                                    if (~changed_num)
                                        if (min_num==1)
                                            min_num = 2;
                                        else
                                            min_num = 1;
                                        end
                                        
                                        E_step = E_step_fnd;
                                        sbnd_E_srch_act = sbnd_E_srch_fnd;
                                        E_acc_act = E_acc_fnd;
                                        changed_num = 1;
                                    else
                                        subband_E(jj+1) = E_result;
                                        %disp('Solution not found, assumed calculated value');
                                        
                                        err = E_result - min_E - sbnd_E_srch_act;
                                        jj = jj + 1;
                                        
                                        sbnd_E_srch_old = sbnd_E_srch_act;
                                        sbnd_E_srch_act = min([max(sbnd_E_srch_act, abs(5*err)), subband_E_search]);
                                        
                                        if (sbnd_E_srch_act<sbnd_E_lim)
                                            sbnd_E_srch_act = sbnd_E_lim;
                                        end
                                        
                                        sbnd_E_lim = max([sbnd_E_lim/10, abs(5*err), min(E_acc,E_acc_act*1e4)]);
                                        
                                        E_step_old = E_step;
                                        E_step = sbnd_E_srch_act/10;
                                        
                                        E_acc_prev = E_acc_act;
                                        E_acc_act = E_acc_act*1e2;
                                        
                                        [f_l, f_h] = PlotWavefunction(E_result, k_t, 1, a);
                                        
                                        %                             figure;
                                        %                             plot(f_h); hold; plot(f_l, 'r'); drawnow;
                                        
                                        accuracy_modifier = min(accuracy_modifier*2, 100);
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Saving the results
                E_k_v{nr} = subband_E;
                k_t_len = [k_t_len, length(k_t_mat)];
                rtsTE(:,:,1) = rtsTE(:,:,2);
                rtsTM(:,:,1) = rtsTM(:,:,2);
                
                t(nr) = toc;
                %disp(['Elapsed time = ', num2str(t(nr)), ' sec']);
                
                % Plotting
                figure(h_E_k); hold on; box on;
                plot(k_t_mat, subband_E/meV);
                title('E-k');
                xlabel('k_t (m^-^1)'); ylabel('E (meV)');
                drawnow;
                
                figure(h_rts_TE); hold on; box on;
                plot(k_t_mat, squeeze(rtsTE(nr,1,1:size(k_t_mat,2))),cl(1));
                for (cb_num=1:num_cb_subbands)
                    plot(k_t_mat, squeeze(rtsTE(nr,cb_num,1:size(k_t_mat,2))),cl(cb_num));
                end
                title('Transmition Strength - k_t');
                ylabel('Transmition Strength'); xlabel('k_t (m^-^1)');
                
                figure(h_rts_TM); hold on; box on;
                plot(k_t_mat, squeeze(rtsTM(nr,1,1:size(k_t_mat,2))),cl(1));
                for (cb_num=1:num_cb_subbands)
                    plot(k_t_mat, squeeze(rtsTM(nr,cb_num,1:size(k_t_mat,2))),cl(cb_num));
                end
                title('Transmition Strength - k_t');
                ylabel('Transmition Strength'); xlabel('k_t (m^-^1)');
                
            end
            
            % Conduction band parabolic dispersion curves
            for (cb_num = 1:num_cb_subbands)
                E_k_c{cb_num} = Materials{2}.E_g*Consts.e_0 + E_0_c(cb_num).*1e-3*Consts.e_0 + (Consts.hbar^2*k_t_mat.^2)./(2*Materials{2}.m_e);
            end
            
        case '4X4_k.p'
            
            for (nk=1:length(k_t_mat))
                %k(nk)=(nk-1)/500; % in A^-1
                %                 
                %l=0; m=1;
                %lm=sqrt((l^2)+(m^2));
                %kx=(l/lm)*k(nk)*1e10;
                %ky=(m/lm)*k(nk)*1e10;
                %k2=(kx^2)+(ky^2);
                %k_t = sqrt(k2);
                %k_t_vec(nk) = k_t;
                k_t = k_t_mat(nk);
                
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
                        H_diag  = [C(zz)*k_t^2+(2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2               ;
                                            sqrt(3)*G2(zz)*k_t^2               , A(zz)*k_t^2+(2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];
                        
                        H_off_p = [       -(D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz))/(2*dz)  ;
                                    sqrt(3)*k_t*(G3(zz))/(2*dz)   ,            -(B(zz))/(2*dzz)        ];
                        
                        H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)        ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                                   -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz) ,            -(B(zz)+B(zz-1))/(2*dzz)        ];
                        
                        H(2*Np-1:2*Np,2*Np-1:2*Np) = H_diag;
                    else
                        H_diag  = [C(zz)*k_t^2+(D(zz+1)+2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                       ;
                                            sqrt(3)*G2(zz)*k_t^2                       , A(zz)*k_t^2+(B(zz+1)+2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];
                        
                        H_off_p = [       -(D(zz+1)+D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ;
                                    sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)        ];
                        
                        H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)         ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                                   -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ,            -(B(zz)+B(zz-1))/(2*dzz)        ];
                        
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
                    F_U_norm = sqrt(trapz(z_grid, F1{nk}(:,ii).^2+F2{nk}(:,ii).^2));
                    Fhh{nk}(:,ii) = F1{nk}(:,ii)./F_U_norm;
                    Flh{nk}(:,ii) = F2{nk}(:,ii)./F_U_norm;
                end
                
                % Overlap integrals --------------------------------------------
                for (cb_index=1:num_cb_subbands)
                    for (vb_index=1:length(E_v{1}))
                        overlap1(nk,cb_index,vb_index) = abs(trapz(z_grid, conj(wf_e_mat{cb_index}(:,2)).*Fhh{nk}(:,vb_index))).^2;
                        overlap2(nk,cb_index,vb_index) = abs(trapz(z_grid, conj(wf_e_mat{cb_index}(:,2)).*Flh{nk}(:,vb_index))).^2;
                        
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
                    hold on; h_F1(ii) = plot(z_grid/1e-10, Fhh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
                    l_F1{ii} = sprintf('F_h_h - h%d',ii);  ylabel('F_h_h');
                    ylabel('F_h_h');
                    subplot(212); grid on; box on;
                    hold on; h_F2(ii) = plot(z_grid/1e-10, Flh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
                    l_F2{ii} = sprintf('F_l_h - h%d',ii); ylabel('F_l_h');
                    ylabel('F_l_h');
                end
                subplot(211);
                hold on; plot(z_grid/1e-10, wf_e_mat{cb_index}(:,2), 'r--', 'linewidth', 2); hold off;
                legend(h_F1,l_F1); xlabel('z (A)');
                title(['Envelope function amplitudes, k_t=',num2str(k_t/1e10),' A^-^1']);
                subplot(212);
                hold on; plot(z_grid/1e-10, wf_e_mat{cb_index}(:,2), 'r--', 'linewidth', 2); hold off;
                legend(h_F2,l_F2); xlabel('z (A)');
                
            end
            
        otherwise
            disp('Unrecognized calculation method');
    end
    
end

%% 7. Saving Results

if (run_cell == 6 || run_cell == -1)
    
    for (ii=1:length(E_k_v))
        E_k_v{ii} = E_k_v{ii}(1:min(k_t_len));
    end
    rtsTE = rtsTE(:,:,1:min(k_t_len));
    rtsTM = rtsTM(:,:,1:min(k_t_len));
    
    Results.N = N;
    Results.Materials = Materials;
    Results.z_grid = z_grid;
    Results.x_profile = x_profile;
    Results.v_e_profile = v_e_profile;
    Results.v_h_profile = v_h_profile;
    Results.m_e_profile = m_e_profile;
    Results.m_h_profile = m_h_profile;
    Results.E_index = E_index;
    Results.wf_e_mat = wf_e_mat;
    Results.wf_h_mat = wf_h_mat;
    Results.k_t_mat = k_t_mat;
    Results.E_0_v = E_state;
    Results.E_k_v = E_k_v;
    Results.E_0_c = E_0_c;
    Results.E_k_c = E_k_c;
    Results.rtsTE = rtsTE;
    Results.rtsTM = rtsTM;
    
end

