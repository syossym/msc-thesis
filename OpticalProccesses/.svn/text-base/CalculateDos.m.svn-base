function OutParams = CalculateDos(InParams)

%
% This function calculates the DOS the quantum structure conduction and valence bands.
%
%   Input: 'InParams' - structure with the input parameters for the calculation.
%
%   Output: 'OutParams' - structure containing the major results of the
%                         simulation.
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, October 2009
% Edited by: -
%

run_cell = -1;    % marks the number of the cell to run.
% Choose '-1' for the execution of the entire script.

%% 1. Constants and definitions

if (run_cell == 1 || run_cell == -1)
    
    global Consts;
    
    N = InParams.N;
    Materials = InParams.Materials
    z_grid = InParams.z_grid;
    m_eff_e = InParams.m_e_profile(length(InParams.m_e_profile)/2);
    V_vb = InParams.v_e_profile;
    v_h_profile = InParams.v_h_profile;
    E_v = InParams.E_k_v;
    E_c = InParams.E_0_c;
    rtsTE = InParams.rtsTE;
    rtsTM = InParams.rtsTM;
    k_t_vec = InParams.k_t_mat(1:length(E_v{1}));
    [vb_len,vc_len] = size(rtsTE(:,:,1));
    
    meV = 1e-3*Consts.e_0;
    
    % Creating the plots
    h_dos_c = figure('Name', 'Conduction Band DOS');
    h_dos_v = figure('Name', 'Valence Band DOS'); 
end

%% 2. Calculation

if (run_cell == 2 || run_cell == -1)
    
    % Conduction band -----------------------------
    
    copyfile('.\Bin\Output\v_e.r', '.\Bin\v.r');
    copyfile('.\Bin\Output\Ee.r', '.\Bin\Ee.r');
    delete('.\Bin\rho.r');
    RunExternalCode('.\Bin\dos.exe', ' -p e');
    copyfile('.\Bin\rho.r', '.\Bin\Output\rho_e.r');
    
    rho_file = load('.\Bin\Output\rho_e.r');
    Dos_e_E_grid = rho_file(:,1);
    Dos_e_bulk = rho_file(:,2);
    Dos_e_qw = rho_file(:,3);
    
    % Valence band --------------------------------
    
    E_grid_step = 0.1e-3;                                % eV
    k_t_step = k_t_vec(2)-k_t_vec(1);
    
    E_v_k = zeros(length(E_v{1}), length(k_t_vec));
    for (vb_index=1:length(E_v))
        for (nk = 1:length(k_t_vec))
            E_v_k(vb_index, nk) = abs(E_v{vb_index}(nk))./Consts.e_0;   % eV
        end
    end
    
    E_min = min(abs(V_vb))./Consts.e_0;                  % eV
    E_max = max(abs(V_vb))./Consts.e_0;                  % eV
    E_grid_DOS = [E_min:E_grid_step:E_max];              % energy grid for the DOS (eV)
    E_grid_size = size(E_grid_DOS,2);
    Dos_E = zeros(length(E_v), length(E_grid_DOS));
    dE_k = diff(E_v_k,1,2);                              % eV
    dk_t = diff(k_t_vec);                                % m
    k_t_mid = k_t_vec(2:length(k_t_vec)) - k_t_step/2;
    E_v_k_node_index = (E_v_k-E_min)./E_grid_step + 1;
    
    for (vb_index = 1:length(E_v))
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

%% 3. Plotting

if (run_cell == 3 || run_cell == -1)
    
    % Conduction band
    figure(h_dos_c); box on; 
    subplot(211); plot(Dos_e_E_grid, Dos_e_bulk);
    title('Conduction band density of states (Dos)');
    ylabel('3D Dos [m^{-2}J^{-1}]'); grid on;
    subplot(212); plot(Dos_e_E_grid, Dos_e_qw, 'g'); 
    xlabel('E [meV]'); ylabel('2D Dos [m^{-2}J^{-1}]'); 
    grid on;
    
    % Valence band
    figure(h_dos_v); box on; grid on;
    Dos_t = zeros(1,length(Dos_E(1,:)));
    for (vb_index=1:vb_len)
        hold on;
        h_dos_fig(vb_index) = plot(E_grid_DOS*1e3, Dos_E(vb_index,:));%, 'Color', colors(vb_index,:));
        Dos_t = Dos_t + Dos_E(vb_index,:);
        hold off;
        l_dos{vb_index} = sprintf('h%d',vb_index);
    end
%     hold on;
%     plot(E_grid_DOS*1e3, Dos_t,':','linewidth',3);
%     hold off;
    xlabel('E [meV]'); ylabel('Dos [m^{-2}eV^{-1}]');
    title('Valence band density of states (DOS)');
    legend(h_dos_fig,l_dos);
    
end

%% 4. Saving Results

if (run_cell == 4 || run_cell == -1)
    
    OutParams.Dos_e_E_grid = Dos_e_E_grid/1e3;   % [eV]
    OutParams.Dos_e_bulk = Dos_e_bulk; 
    OutParams.Dos_e_qw = Dos_e_qw;               % [eV^-1 m^-2)]
    OutParams.Dos_h_E_grid = E_grid_DOS;         % [eV]
    OutParams.Dos_h_k = Dos_k;
    OutParams.Dos_h_E = Dos_E;
   
end