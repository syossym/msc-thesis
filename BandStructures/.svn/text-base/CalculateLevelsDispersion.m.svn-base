function [Cond, Valence, QStruct, OpticalParams] = CalculateLevelsDispersion(QStruct, Cond, Valence, params)
%
% This function calculates the dispersion relations of the subbands in the
% conduction and valence bands.
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'Cond' - calculated conduction band parameters.
%          'Valence' - calculated valence band parameters.
%          'params' - simulaion parameters.
%
%   Output: 'Cond' - calculated conduction band parameters.
%           'Valence' - calculated valence band parameters.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Defintions

k_t_max = 1e9;                    % max. calculated k_t
k_t_step = 1e7;                   % k_t step for E-k dispersion
k_t_vec = [0:k_t_step:k_t_max];
params.k_t_vec = k_t_vec;
QStruct.k_t_vec = k_t_vec;

%% Conduction band

for (cc=1:params.num_cond_subbands)
    Cond.E_k(cc,:) = Cond.E_0(cc)*Consts.e_0 + (Consts.hbar^2*k_t_vec.^2)./(2*QStruct.ActiveLayers{1}.m_e);
end

%% Valence band

switch (params.ValenceCalculationMethod)
    case 'PropagationMatrix',
        [Cond, Valence, QStruct, OpticalParams] = PropagatingMatrixMethod(QStruct, Valence, Cond, params);
    case 'k.p',
        [Valence, OpticalParams] = KdotPMethod(QStruct, Valence, Cond, params);
end

function Profile = DecimatePotentialProfile(L_clad, QStruct)

global Consts;

dz = QStruct.z_grid(2)-QStruct.z_grid(1);

L = 0;
for (ll=1:length(QStruct.Layers))
    if (QStruct.Layers{ll}.IsActive)
        z = L-L_clad:dz:L+QStruct.Layers{ll}.L*1e-10+L_clad;
        z = z(1:end-1);
        Profile.v_h_profile = interp1(QStruct.z_grid, -(QStruct.Vh_profile-QStruct.Phi), z, 'cubic').'*Consts.e_0;
        Profile.v_h_profile = Profile.v_h_profile + abs(min(Profile.v_h_profile));
        Profile.gamma1_profile = interp1(QStruct.z_grid, QStruct.gamma1_profile, z, 'cubic').';
        Profile.gamma2_profile = interp1(QStruct.z_grid, QStruct.gamma2_profile, z, 'cubic').';
        Profile.gamma3_profile = interp1(QStruct.z_grid, QStruct.gamma3_profile, z, 'cubic').';
        Profile.z_grid_partial = z;
        Profile.z_grid = (0:dz:dz*(length(Profile.z_grid_partial)-1)).';
        break;
    end
    L = L + QStruct.Layers{ll}.L*1e-10;
end

function [Cond, Valence, QStruct, OpticalParams] = PropagatingMatrixMethod(QStruct, Valence, Cond, params)

global Consts;

%% Simulation defintions

L_clad = 200e-10;    % [m]

Decimated_Profile = DecimatePotentialProfile(L_clad, QStruct);
QStruct.ActiveLayers{1}.z_grid = Decimated_Profile.z_grid;

k_t_mat = params.k_t_vec;
k_t_step = k_t_mat(2)-k_t_mat(1);
k_t_len = [];

dz = Decimated_Profile.z_grid(2)-Decimated_Profile.z_grid(1);      % grid step size [m]
N  = size(Decimated_Profile.v_h_profile,1);                        % num of grid points
nrm = 1e7;                                                         % normalization coefficient
meV = 1e-3*Consts.e_0;                                             % [J]

min_E_diff = 0.05*meV;         % helps to determine whether two solutions are the same [J]
a_lim = 1e5;                   % max. for a
E_min = min(Decimated_Profile.v_h_profile);      % energy minimum [J]
E_max = max(Decimated_Profile.v_h_profile); % min(v_h_profile(1),...
%v_h_profile(N));   % energy maximum [J]
E_step = 0.1*meV;    	       % energy step for band-edge calculation [J]
E_acc  = 1e-9*meV;             % energy resolution of the calculated subband energies [J]

% Search parameters
subband_E_search  = 0.15*meV;  % energy search range [J]
E_abort           = meV;       % abort search when E range > E_abort

%% Band-edge calculation

h_LH_Subbnads = figure('Name','LH Subbnads');
h_HH_Subbnads = figure('Name','HH Subbnads');
h_E_k = figure('Name','E_k');
h_rts_TE = figure('Name','RTS TE');
h_rts_TM = figure('Name','RTS TM');
cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';...
    'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];
ch = cl;

% Prepare the plots
figure(h_LH_Subbnads);
plot(Decimated_Profile.z_grid/1e-10, Decimated_Profile.v_h_profile/meV, 'k');
hold; title('LH Subbands'); xlabel('z [A]'); ylabel('E [meV]');
figure(h_HH_Subbnads);
plot(Decimated_Profile.z_grid/1e-10, Decimated_Profile.v_h_profile/meV, 'k');
hold; title('HH Subbands'); xlabel('z [A]'); ylabel('E [meV]');

% Start energy loop
nr = 1;
old_sign_l = 1;
old_sign_h = 1;
old_sign_d = 1;
old_sign_d2 = 1;
old_sign_d3 = 1;
ii = 1;
k_t = 1e5;

for (E = E_min:E_step:E_max)
    %disp(['Step ' num2str(ii) ': E=' num2str(E)]);

    % Transfer matrix
    tf_mat = TransferMatrix(E, k_t, Decimated_Profile);
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

                [a, E_result] = WavefunctionMinimum(Decimated_Profile,E-2*E_step,E,k_t,E_acc);

                %disp(['E_result=' num2str(E_result/meV)]);

                if (nr==1) | ( (nr>1) & abs(E_result-E_state(nr-1))>min_E_diff )
                    E_state(nr) = E_result;
                    [f_h, f_l]  = PlotWavefunction(E_result,k_t,1,a,Decimated_Profile);

                    %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);

                    if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4 )
                        E_state(nr) = E_result;
                        coeff(nr)   = a;
                        wf_l(nr,:)  = f_l;
                        wf_h(nr,:)  = f_h;

                        % Plotting
                        figure(h_LH_Subbnads);
                        plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                        drawnow;
                        refresh(h_LH_Subbnads);
                        figure(h_HH_Subbnads);
                        plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                        drawnow;
                        refresh(h_HH_Subbnads);

                        wf_h_mat{nr} = [Decimated_Profile.z_grid, f_h', f_l'];

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

                [a, E_result] = WavefunctionMinimum(Decimated_Profile,E-2*E_step,E+2*E_step,k_t,E_acc,E_step/50);

                if (nr==1) | ( (nr>1) & abs(E_result-E_state(nr-1))>min_E_diff )

                    E_state(nr) = E_result;
                    [f_h,f_l]   = PlotWavefunction(E_result,k_t,1,a,Decimated_Profile);

                    %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);

                    if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4  )
                        E_state(nr) = E_result;
                        coeff(nr)   = a;
                        wf_l(nr,:)  = Normalize(f_l);
                        wf_h(nr,:)  = Normalize(f_h);

                        % Plotting
                        figure(h_LH_Subbnads);
                        plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                        drawnow;
                        refresh(h_LH_Subbnads);
                        figure(h_HH_Subbnads);
                        plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                        drawnow;
                        refresh(h_HH_Subbnads);

                        wf_h_mat{nr} = [Decimated_Profile.z_grid, f_h', f_l'];

                        nr = nr + 1;

                    else

                        [a, E_result] = WavefunctionMinimum(Decimated_Profile,E_result,E+2*E_step,k_t,E_acc,E_step/50,2);
                        E_state(nr) = E_result;
                        [f_h,f_l]   =  PlotWavefunction(E_result,k_t,1,a,Decimated_Profile);

                        %disp(['Edge=' num2str((f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2))]);

                        if ( (f_l(N)^2 + f_h(N)^2)/(max(abs(f_l))^2+max(abs(f_h))^2) < 1e-4 )
                            E_state(nr) = E_result;
                            wf_l(nr,:)  = Normalize(f_l);
                            wf_h(nr,:)  = Normalize(f_h);

                            % Plotting
                            figure(h_LH_Subbnads);
                            plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_l(nr,:).^2)/nrm, cl(nr,:));
                            drawnow;
                            refresh(h_LH_Subbnads);
                            figure(h_HH_Subbnads);
                            plot(Decimated_Profile.z_grid/1e-10, E_state(nr)/meV + (wf_h(nr,:).^2)/nrm, ch(nr,:));
                            drawnow;
                            refresh(h_HH_Subbnads);

                            wf_h_mat{nr} = [Decimated_Profile.z_grid, f_h', f_l'];

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

% Calculate the optical transition strength
for (vb_num=1:params.num_valence_subbands)
    Valence.Wf_0_Q{vb_num} = wf_h_mat{vb_num};
    for (cb_num=1:params.num_cond_subbands)
        wf_c =  Normalise_Wfs(Decimated_Profile.z_grid, interp1(QStruct.z_grid, Cond.Wf_0(cb_num,:), Decimated_Profile.z_grid_partial, 'pchip'), 'Rel');
        f_h_norm = wf_h_mat{vb_num}(:,2).';
        f_l_norm = wf_h_mat{vb_num}(:,3).';
        overlap1 = trapz(Decimated_Profile.z_grid, (wf_c.*conj(f_h_norm)));
        overlap2 = trapz(Decimated_Profile.z_grid, (wf_c.*conj(f_l_norm)));
        rtsTE(vb_num, cb_num, 1) = (3/4)*(abs(overlap1^2) + (1/3)*abs(overlap2^2));
        rtsTM(vb_num, cb_num, 1) = abs(overlap2^2);
    end
end


%% Dispersion calculation

for (vv=1:params.num_valence_subbands)

    % Definitions
    E_0 = E_state(vv);    % subband edge energy [J]
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
    
    tic;

    % Iterating over k_t
    while ((jj < size(k_t_mat,2)) && ~abort)

        k_t = k_t_mat(jj+1);

        min_E = subband_E(jj) + dE - sbnd_E_srch_act;
        max_E = subband_E(jj) + dE + sbnd_E_srch_act;
        min_E_arc(jj+1) = min_E;
        max_E_arc(jj+1) = max_E;

        [a,E_result,err] = WavefunctionMinimum(Decimated_Profile,min_E,max_E,k_t,E_acc_act,E_step,min_num);

        if (jj/50==floor(jj/50))
            drawnow;
        end

        if (err)    % -- minimum not found

            if ((min_num==2) && (jj<4))
                min_num = 1;
            else
                num_steps = 2*sbnd_E_srch_act/E_step;

                if ((E_step>E_step_old/10) && (sbnd_E_srch_act<sbnd_E_srch_old*10))

                    if (num_steps<3e3)
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

            [f_h,f_l] = PlotWavefunction(E_result, k_t, 1, a, Decimated_Profile);

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

                %                 figure(1);
                for (cb_num=1:Cond.num_subbands)
                    wf_c =  Normalise_Wfs(Decimated_Profile.z_grid, interp1(QStruct.z_grid, Cond.Wf_0(cb_num,:), Decimated_Profile.z_grid_partial, 'pchip'), 'Rel');
                    Cond.Wf_0_Q(cb_num,:) = wf_c;
                    %f_h_norm = Normalise_Wfs(Decimated_Profile.z_grid,f_h,'Abs');
                    %f_l_norm = Normalise_Wfs(Decimated_Profile.z_grid,f_l,'Abs');
                    f_h_norm = f_h;
                    f_l_norm = f_l;
                    overlap1 = trapz(Decimated_Profile.z_grid, (wf_c.*conj(f_h_norm)));
                    overlap2 = trapz(Decimated_Profile.z_grid, (wf_c.*conj(f_l_norm)));

                    %                     subplot(211); plot(Decimated_Profile.z_grid, wf_c, 'r', Decimated_Profile.z_grid, f_h_norm, 'b');
                    %                     subplot(212); plot(Decimated_Profile.z_grid, wf_c, 'r', Decimated_Profile.z_grid, f_l_norm, 'b');
                    %
                    rtsTE(vv, cb_num, jj+1) = (3/4)*(abs(overlap1^2) + (1/3)*abs(overlap2^2));
                    rtsTM(vv, cb_num, jj+1) = abs(overlap2^2);

                    %rtsTE(vv, cb_num, jj+1) = rtsTE(vv, cb_num, jj+1)/3;
                    %rtsTM(vv, cb_num, jj+1) = rtsTM(vv, cb_num, jj+1)/3;

                    %rtsTE(vv, cb_num, jj+1) = 0.5*(abs(overlap1^2) + (1/3)*abs(overlap2^2));
                    %rtsTM(vv, cb_num, jj+1) = (2/3)*abs(overlap2^2);
                end

                err = E_result - min_E - sbnd_E_srch_act;
                jj = jj + 1;

                %disp(['k_t=', num2str(k_t), ', E_acc=', num2str(E_acc_act/meV), ...
                %', srch_dE=', num2str(sbnd_E_srch_act/meV), ' E_step=', num2str(E_step/meV) ...
                %' E found ', num2str(E_result/meV), ' meV']);

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
    E_k_v{vv} = subband_E;
    Valence.E_k(vv,:) = interp1(k_t_mat, subband_E, params.k_t_vec, 'pchip');
    for (cc=1:params.num_cond_subbands)
        OpticalParams.rts_TE(vv, cc, :) = interp1(k_t_mat, squeeze(rtsTE(vv,cc,:)), params.k_t_vec, 'pchip');
        OpticalParams.rts_TM(vv, cc, :) = interp1(k_t_mat, squeeze(rtsTM(vv,cc,:)), params.k_t_vec, 'pchip');
    end
    k_t_len = [k_t_len, length(k_t_mat)];

    t(vv) = toc;
    %disp(['Elapsed time = ', num2str(t(vv)), ' sec']);

    % Plotting
    figure(h_E_k); hold on; box on;
    plot(k_t_mat, subband_E/meV);
    title('E-k');
    xlabel('k_t (m^-^1)'); ylabel('E (meV)');
    drawnow;

    figure(h_rts_TE); hold on; box on;
    %plot(k_t_mat, squeeze(rtsTE(vv,1,1:size(k_t_mat,2))),cl(1));
    for (cb_num=1:Cond.num_subbands)
        plot(k_t_mat, (squeeze(rtsTE(vv,cb_num,1:size(k_t_mat,2)))),cl(cb_num));
    end
    title('TE');
    ylabel('Transmition Strength'); xlabel('k_t (m^-^1)');

    figure(h_rts_TM); hold on; box on;
    %plot(k_t_mat, squeeze(rtsTM(vv,1,1:size(k_t_mat,2))),cl(1));
    for (cb_num=1:Cond.num_subbands)
        plot(k_t_mat, (squeeze(rtsTM(vv,cb_num,1:size(k_t_mat,2)))),cl(cb_num));
    end
    title('TM');
    ylabel('Transmition Strength'); xlabel('k_t (m^-^1)');

end

function [Valence, OpticalParams] = KdotPMethod(QStruct, Valence, Cond, params)

global Consts;

OpticalParams = [];

%% Simulation parameters

%L_clad = 100e-10;    % [m]
%Profile = DecimatePotentialProfile(L_clad, QStruct);
Profile = QStruct;

G1 = Profile.gamma1_profile;
G2 = Profile.gamma2_profile;
G3 = Profile.gamma3_profile;
A = G1-G2;
B = G1+2*G2;
C = G1+G2;
D = G1-2*G2;
E = G2+G3;
V = -(Profile.Vh_profile-Profile.Phi).*Consts.e_0;

Np = length(Profile.z_grid);
k_t_vec = params.k_t_vec;
dz = Profile.z_grid(2)-Profile.z_grid(1);
dzz = dz^2;

%% k.p calculation

%Valence.E_k(1,:) = Valence.E_0(1:5);

switch (params.KdotPType)
    case '4x4',

        for nk=1:length(k_t_vec)
            %             k(nk)=(nk-1)/500; % in A^-1
            %             l=0;m=1;lm=sqrt((l^2)+(m^2));
            %             kx=(l/lm)*k(nk)*1e10;ky=(m/lm)*k(nk)*1e10;
            %             k2=(kx^2)+(ky^2);
            %             k_t = sqrt(k2);
            %             k_t_vec(nk) = k_t;
            k_t = k_t_vec(nk)

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
            %[nk sum(sum(abs(H-H')))];
            H_sparse = sparse(H);
            clear H;
            opt.disp = 0;
            [F,En] = eigs(H_sparse, 5, 0, opt);
            %[F,En] = eig(H);
            En = diag(En)./Consts.e_0;
            [En,I] = sort(real(En));
            En = -En;
            F1{nk} = F(1:2:end,I);
            F2{nk} = F(2:2:end,I);

            E_v{nk} = En(1:5);

            E1(nk)=En(1); E2(nk)=En(2); E3(nk)=En(3); E4(nk)=En(4);
            E5(nk)=En(5);

            Valence.E_k(nk,:) = En(1:5);

            %             for (ii=1:length(E_v{1}))
            %                 F_U_norm = sqrt(trapz(z, F1{nk}(:,ii).^2+F2{nk}(:,ii).^2));
            %                 Fhh{nk}(:,ii) = F1{nk}(:,ii)./F_U_norm;
            %                 Flh{nk}(:,ii) = F2{nk}(:,ii)./F_U_norm;
            %             end
            %
            %             % Overlap integrals --------------------------------------------
            %             for (cb_index=1:length(Cond.E_0))
            %                 for (vb_index=1:length(E_v{1}))
            %                     overlap1(nk,cb_index,vb_index) = trapz(z, Cond.Wf_0(cb_index,:)'.*Fhh{nk}(:,vb_index)).^2;
            %                     overlap2(nk,cb_index,vb_index) = trapz(z, Cond.Wf_0(cb_index,:)'.*Flh{nk}(:,vb_index)).^2;
            %
            %                     rtsTE(nk,cb_index,vb_index) = (3/4)*overlap1(nk, cb_index,vb_index) + ...
            %                                                   (3/12)*overlap2(nk, cb_index,vb_index);
            %                     rtsTE_k(vb_index,cb_index,nk) = rtsTE(nk,cb_index,vb_index);
            %                     rtsTM(nk,cb_index,vb_index) = overlap2(nk, cb_index,vb_index);
            %                     rtsTM_k(vb_index,cb_index,nk) = rtsTM(nk,cb_index,vb_index);
            %                 end
            %             end
            %
            %             % Plot the envelope functions ----------------------------------
            %             figure(h_env_func);
            %             clf;
            %             colors = [1 0 0; 1 1 1; 0 1 0; 0 1 1];
            %             for(ii=1:length(E_v{1}))
            %                 subplot(211); grid on; box on;
            %                 hold on; h_F1(ii) = plot(z/1e-10, Fhh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
            %                 l_F1{ii} = sprintf('F_h_h - h%d',ii); ylabel('F_h_h');
            %                 ylabel('F_h_h');
            %                 subplot(212); grid on; box on;
            %                 hold on; h_F2(ii) = plot(z/1e-10, Flh{nk}(:,ii), 'Color', colors(ii,:)); hold off;
            %                 l_F2{ii} = sprintf('F_l_h - h%d',ii); ylabel('F_l_h');
            %                 ylabel('F_l_h');
            %             end
            %             subplot(211);
            %             hold on; plot(z/1e-10, Cond.Wf_0, 'r--', 'linewidth', 2); hold off;
            %             legend(h_F1,l_F1); xlabel('z (A)');
            %             title(['Envelope function amplitudes, k_t=',num2str(k_t/1e10),' A^-^1']);
            %             subplot(212);
            %             hold on; plot(z/1e-10, Cond.Wf_0, 'r--', 'linewidth', 2); hold off;
            %             legend(h_F2,l_F2); xlabel('z (A)');

        end
end



