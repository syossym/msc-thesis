function Results = Bulk_kp(mat_name, grid_size)

%
% This function calculates the band structure according to the k.p for
% a bulk semiconductor.
%
%   Input: 'mat_name'  - the name of the SC material
%          'grid_size' - the size of the simualation grid
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

    % Local grid definitions
    %z = zeros(5);
    %Z = zeros(10);

    % Defining the material parameters
    Params.T = 300; Params.E = 0;
    Mat = GetMaterial(mat_name, Params);

    t1 = (Consts.hbar^2)*Mat.g1/(2*Consts.m_0*Consts.e_0*(Mat.a^2));
    t2 = (Consts.hbar^2)*Mat.g2/(2*Consts.m_0*Consts.e_0*(Mat.a^2));
    t3 = (Consts.hbar^2)*Mat.g3/(2*Consts.m_0*Consts.e_0*(Mat.a^2));

    Nt = grid_size;
    kk = 1*linspace(0,1,Nt);   % k*a
    Ev = 0;
    
    % Crystalographic orientations, [n,m,l]
    lattice_vec = zeros(2,3);        
    lattice_vec(1,:) = [1, 0, 0];        % X-direction 
    lattice_vec(2,:) = [1, 1, 1];        % S-direction  

    k_vec = zeros(length(lattice_vec(:,1)),length(kk));
    
end

%% 2. Simulation

if (run_cell == 2 || run_cell == -1)

    for (ll = 1:length(lattice_vec(:,1)))
        for (Nk = 1:Nt)
            k = 2*pi*kk(Nk)*lattice_vec(ll,:);
            k_vec(ll,Nk) = norm(k);

            % Luttinger-Kohn model parameters
            P  = Ev + (t1*sum(k.*k));
            Q  = t2*((k(1)^2) + (k(2)^2) - (2*(k(3)^2)));
            R  = -(sqrt(3)*t2*((k(1)^2)-(k(2)^2)))+(i*2*t3*sqrt(3)*k(1)*k(2));
            S  = 2*t3*sqrt(3)*((k(1)-(i*k(2)))*k(3));

            % 4X4
            H4 = -[P+Q -S    R   0;
                   -S' P-Q   0   R;
                    R'  0   P-Q  S;
                    0   R'   S' P+Q];

            [V,D]  = eig(H4);
            eiglst = sum(D);
            ELK4(Nk,:) = sort(real(eiglst));
            
            % 6X6
            H6  =  -[P+Q        -S            R           0       -S/sqrt(2)       sqrt(2)*R;
                     -S'        P-Q           0           R       -sqrt(2)*Q       sqrt(1.5)*S;
                      R'         0           P-Q          S        sqrt(1.5)*S'    sqrt(2)*Q;
                      0          R'           S'         P+Q      -sqrt(2)*R'     -S'/sqrt(2);
                  -S'/sqrt(2) -sqrt(2)*Q'  sqrt(1.5)*S -sqrt(2)*R  P+Mat.Del          0;
                   sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q'  -S/sqrt(2)       0          P+Mat.Del];

            [V,D]  = eig(H6);
            eiglst = sum(D);
            ELK6(Nk,:) = sort(real(eiglst));

            % 8X8
            Pz = i*(P/sqrt(3));
            Pp = Pz*(k(1)+k(2));
            Pm = Pz*(k(1)+k(2));
            param = Mat.E_g + 8*(Consts.hbar^2/(2*Consts.m_0))*norm(k)^2;
            
            H8  =  [     param        -sqrt(2)*Pz*k(3)     Pz*k(3)     sqrt(3/2)*Pp      0      -Pm/sqrt(2)     -Pm             0     ;
                     -sqrt(2)*Pz'*k(3)     -(P-Q)         -sqrt(2)*Q         S'      -Pp/sqrt(2)       0      -sqrt(3/2)*S       R     ;
                         Pz'*k(3)         -sqrt(2)*Q      -(P+Mat.Del)  -S'/sqrt(2)     -Pp'    sqrt(3/2)*S      0          sqrt(2)*R ;
                       sqrt(3/2)*Pp'          S            -S/sqrt(2)     -(P+Q)          0           -R       -sqrt(2)*R        0     ;
                           0             -Pp/sqrt(2)          -Pp            0          param   sqrt(2)*Pz*k(3)  -Pz*k(3)    -sqrt(3/2)*Pm ;
                       -Pm'/sqrt(2)           0            sqrt(3/2)*S'     -R'     sqrt(2)*Pz'*k(3) -(P-Q)     -sqrt(2)*Q       S     ; 
                         -Pm'           -sqrt(3/2)*S'          0         -sqrt(2)*R'   -Pz'*k(3)   -sqrt(2)*Q  -(P-Mat.Del)   -S/sqrt(2) ;
                           0                  R'            sqrt(2)*R'       0       -sqrt(3/2)*Pm'     S'      -S'/sqrt(2)    -(P+Q)    ];

            [V,D]  = eig(H8);
            eiglst = sum(D);
            ELK8(Nk,:) = sort(real(eiglst));
            
        end
        
        Eigs.kp4{ll} = ELK4;
        Eigs.kp6{ll} = ELK6;
        Eigs.kp8{ll} = ELK8;
    end

end

%% 3. Plotting

if (run_cell == 3 || run_cell == -1)

    figure(1); hold on;
    figure(2); hold on;

    for (ll = 1:length(lattice_vec(:,1)))
        if (rem(ll,2) == 0)
            k_vec(ll,:) = -k_vec(ll,:);
        end
        
        E_e_k = Mat.E_g*Consts.e_0 + (Consts.hbar^2*(k_vec(ll,:)/Mat.a).^2)./(2*Mat.m_e);
        
        figure(1);
        %plot(k_vec(ll,:)./Mat.a, Eigs.kp4{ll}, 'b');
        plot(k_vec(ll,:)./Mat.a, Eigs.kp6{ll}, 'r-');
        %plot(k_vec(ll,:)./Mat.a, Eigs.kp8{ll}, 'k');
        plot(k_vec(ll,:)./Mat.a, E_e_k/Consts.e_0, 'r');
        xlabel('k (m^{-1})');
        ylabel('Energy (eV)');
        box on;
        axis([-1.2e9, 1.2e9, -0.8, 2]);

        figure(2);
        %plot((k_vec(ll,:)./Mat.a)/(2*pi/Mat.a), Eigs.kp4{ll}, 'b');
        %plot((k_vec(ll,:)./Mat.a)/(2*pi/Mat.a), Eigs.kp6{ll}, 'r-');
        plot((k_vec(ll,:)./Mat.a)/(2*pi/Mat.a), Eigs.kp8{ll}, 'k');
        plot((k_vec(ll,:)./Mat.a)/(2*pi/Mat.a), E_e_k/Consts.e_0, 'r');
        xlabel('k/(2\pi/a)');
        ylabel('Energy (eV)');
        box on;
    end
end

%% 4. Saving Results

if (run_cell == 4 || run_cell == -1)

    Results.k_vec = k_vec;
    Results.Eigs = Eigs;

end

