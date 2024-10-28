function [Cond, Valence] = SolveSchrodingerEquation(QStruct, params)
%
% This function solves the Schrodinger equation and finds the energies
% and wavefunction for the conduction and valence bands of the quantum
% structure.
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'params' - simulaion parameters.
%
%   Output: 'Cond' - structure containig the conduction band energies and
%                    wavefucntion.
%           'Valence' - structure containig the valence band energies and
%                       wavefucntion.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

h_LH_Subbnads = figure('Name','LH Subbnads');
h_HH_Subbnads = figure('Name','HH Subbnads');

%% Conduction band

delete('.\Bin\v.r');
RunExternalCode('.\Bin\efxv.exe', ' -p e');
delete('.\Bin\m.r');
RunExternalCode('.\Bin\efxm.exe', ' -p e');
delete('.\Bin\Ee.r');

% Calculating energy levels
RunExternalCode(['.\Bin\efshoot.exe', ' -p e -s ', num2str(params.cb_subbands)]);
copyfile('.\Bin\Ee.r', '.\Bin\Output\Ee.r');

% Calculating wavefunctions
RunExternalCode('.\Bin\efwf.exe');
copyfile('.\Bin\wf_e*', '.\Bin\Output\');

E_0_c = load('.\Bin\Output\Ee.r');
for (cb_num = 1:num_cb_subbands)
    eval(['wf_e_temp=load(''.\Bin\Output\wf_e' num2str(cb_num) '.r'');']);
    wf_e_mat{cb_num} = wf_e_temp;
end

Cond.num_subbands = params.cb_subbands;
Cond.E_0 = E_0_c;
Cond.Wf_0 = wf_e_mat;

%% Valence band

dz = QStruct.z_grid(2)-QStruct.z_grid(1);      % grid step size [m]
N  = size(QStruct.Vh_profile,1);               % num of grid points
nrm = 1e7;                                     % normalization coefficient
meV = 1e-3*Consts.e_0;                         % [J]

min_E_diff = 0.05*meV;                % helps to determine whether two solutions are the same [J]
a_lim = 1e5;                          % max. for a
k_t = 1e5;                            % k_t for band-edge calculation
E_min = min(QStruct.Vh_profile);      % energy minimum [J]
E_max = max(QStruct.Vh_profile);      % min(v_h_profile(1),...
%v_h_profile(N));   % energy maximum [J]
E_step = 0.1*meV;           	      % energy step for band-edge calculation [J]
E_acc  = 1e-9*meV;                    % energy resolution of the calculated subband energies [J]

% Search parameters
subband_E_search  = 0.15*meV;  % energy search range [J]
E_abort           = meV;       % abort search when E range > E_abort
k_t_max           = 1e9;       % max. calculated k_t
k_t_step          = 1e7;       % k_t step for E-k dispersion
num_cb_subbands   = QStruct.cb_subbands;         % number of conduction band subbands taken into account
%num_vb_subbands   = nr-1;      % number of valence band subbands taken into account

% Plotting parameters
cl = ['m  ';'c  ';'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';...
    'r  ';'g  ';'b  ';'y  ';'m-.';'c-.';'r-.';'g-.';'b-.';'y-.';'m  ';'c  ';'r  ';'g  ';'b  ';'y  '];
ch = cl;

% Prepare the plots
figure(h_LH_Subbnads);
plot(z_grid/1e-10, v_h_profile/meV, 'k');
hold; title('LH Subbands'); xlabel('z [A]'); ylabel('E [meV]');
figure(h_HH_Subbnads);
plot(z_grid/1e-10, QStruct.Vh_profile/meV, 'k');
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

Valence.num_subbands = length(E_state);
Valence.E_0 = E_state;
Valence.Wf_0 = wf_h_mat;
