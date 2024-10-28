clc; close all; warning off;

% Init Simulation
global Consts;
global R G_TE G_TM z_exp;
global project_path;

project_path = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code';
cd(project_path);
run('.\Common\AddPath.m');

% Create the physical constants structure
Constants;
params.x = 0.1;
params.T = 2;
params.E = 0;
M_GaAs = GetMaterial('GaAs', params);
M_AlAs = GetMaterial('AlAs', params);
M_GaAlAs = GetMaterial('GaAlAs', params);

% Definitions
c = 2.99792458e8;  % [m/s]
E_vec = [1.52:0.0001:1.53]*Consts.e_0;
lambda_vec = (2*pi./E_vec)*Consts.c*Consts.hbar; %[1300:0.1:1600].*1e-9;  % [m]
E_c = 1.5244*Consts.e_0;
lambda_c = (2*pi/E_c)*Consts.c*Consts.hbar; %min(lambda_vec)+0.5*(max(lambda_vec)-min(lambda_vec))      % [m]
n_H = M_GaAlAs.n;
n_L = M_AlAs.n;
n_c = n_H;
n_c_l = 1; %M_GaAs.n; %2.98;  % refractive index of the left cladding
n_c_r = 1; %M_GaAs.n;  % refractive index of the right cladding
n_c_c = M_GaAlAs.n;      % cavity QW cladding refractive index
n_c_w = M_GaAs.n;        % cavity QW refractive index
N_l = 15;     % number of alternating index pairs
N_r = 25;
grid = 300;
m = 2;
E_in = [1;0];     % left-hand side incident field amplitude

% Parameter calculation
l_L = lambda_c/n_L/4;
l_H = lambda_c/n_H/4;
l_c = m*lambda_c/n_c/2;
l_c_w = 200e-10;           % QW width [m]
l_c_c = (l_c-l_c_w)/2;     % QW cladding layer widths [m]

n_vec_l = [repmat([n_H, n_L],1,N_l)];
l_vec_l = [repmat([l_H, l_L],1,N_l)];
n_profile_l = [n_c_l.*ones(1,floor(l_L/1e-9)*grid), repmat([n_H*ones(1,floor(l_H/1e-9)*grid), n_L*ones(1,floor(l_L/1e-9)*grid)],1,N_l)];

n_vec_r = [n_L, repmat([n_H, n_L],1,N_r)];
l_vec_r = [l_L, repmat([l_H, l_L],1,N_r)];
n_profile_r = [n_L*ones(1,grid), repmat([n_H*ones(1,floor(l_H/1e-9)*grid), n_L*ones(1,floor(l_L/1e-9)*grid)],1,N_r), n_H*ones(1,floor(l_L/1e-9)*grid)];

n_vec = [n_vec_l, n_c_c, n_c_w, n_c_c, n_vec_r];
n_profile = [n_profile_l, n_c_c*ones(1,floor(l_c_c/1e-9)*grid), n_c_w*ones(1,floor(l_c_w/1e-9)*grid), n_c_c*ones(1,floor(l_c_c/1e-9)*grid), n_profile_r];
l_vec = [l_vec_l, l_c_c, l_c_w, l_c_c, l_vec_r];

dl = sum(l_vec)./length(n_profile);      % [m]
l_profile = [grid*dl*(-1:1/grid:0), (length(n_profile)-grid)*dl*(0:1/(length(n_profile)-grid):1)];

z_grid_calc = 0;
n_profile = n_vec(1);
for (ll=1:length(l_vec))
    vec = z_grid_calc(end)+l_vec(ll)/grid : l_vec(ll)/grid : z_grid_calc(end)+l_vec(ll);
    z_grid_calc = [z_grid_calc, vec ];
    n_profile = [n_profile, ones(1,length(vec))*n_vec(ll)];
end
n_vec_passive = n_vec;

% Simulation
for (con_num=1:length(con_vec))
    dn_w_vec = interp1(E_exc, smooth(delta_n_TE_HF_vec{con_num}), E_vec, 'pchip');
    n_vec = n_vec_passive;
    for (gg=1:length(lambda_vec))
        dn_vec = [zeros(1,length(n_vec_l)), 0, dn_w_vec(gg), 0, zeros(1,length(n_vec_r))];
        n_vec = n_vec_passive + dn_vec;
     
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
            (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
            (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
        
        %         M_DBR = M_c_r*M_DBR;
        %         phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
        %         %phi_vec(gg) = 0;
        %
        %         for (nn=length(n_vec):-1:1)
        %             k = (2*pi/lambda_vec(gg))*n_vec(nn);
        %             %phi_vec(gg) = phi_vec(gg) + (2*pi/lambda_vec(gg))*n_vec(nn)*l_vec(nn);
        %             M_p = [exp(1i*k*l_vec(nn)), 0 ;...
        %                 0 , exp(-1i*k*l_vec(nn))];
        %             if (nn~=1)
        %                 M_i = (1/(2*n_vec(nn))).*[(n_vec(nn)+n_vec(nn-1)) , (n_vec(nn)-n_vec(nn-1)) ; ...
        %                     (n_vec(nn)-n_vec(nn-1)) , (n_vec(nn)+n_vec(nn-1))];
        %             else
        %                 M_i = eye(2);
        %             end
        %
        %             M_DBR = M_i*M_p*M_DBR;
        %
        %         end
        %         M_DBR = M_c_l*M_DBR;
        
        M_DBR = M_DBR*M_c_l;
        phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
        for (nn=1:1:length(n_vec))
            k = (2*pi/lambda_vec(gg))*n_vec(nn);
            %phi_vec(gg) = phi_vec(gg) + (2*pi/lambda_vec(gg))*n_vec(nn)*l_vec(nn);
            M_p = [exp(1i*k*l_vec(nn)),         0          ;...
                0             , exp(-1i*k*l_vec(nn))];
            if (nn~=length(n_vec))
                M_i = (1/(2*n_vec(nn))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
                    (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];
            else
                M_i = eye(2);
            end
            
            M_DBR = M_DBR*M_p*M_i;
        end
        M_DBR = M_DBR*M_c_r;
        r(gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    
    % Plotting
    figure(1);
    subplot(311);
    %plot(l_profile(1:length(n_profile))./1e-9, n_profile); hold on;
    plot(z_grid_calc/1e-10, n_profile); hold on;
    %vec = 0:max(n_profile)/length(n_profile):max(n_profile)-max(n_profile)/length(n_profile);
    %plot(zeros(1,length(vec)), vec, ':g',...
    %    ones(1,length(vec)).*(max(l_profile)-grid*dl)/1e-9, vec, ':g'); hold on;
    %axis([min(l_profile)./1e-9, max(l_profile)./1e-9, 0 , max(n_profile)+0.2]);
    axis([min(z_grid_calc)/1e-10, max(z_grid_calc)/1e-10, 0 , max(n_profile)+0.2]);
    xlabel('z [A]'); ylabel('n');
    title('(a)');
    subplot(312); hold on;
    %plot(lambda_vec./lambda_c, abs(r)); hold on;
    [AX,H1,H2] = plotyy(lambda_vec./lambda_c, abs(r), lambda_vec./lambda_c, angle(r)/pi);
    set(get(AX(2), 'Ylabel'), 'String', '\alpha_r/\pi');
    set(H1,'LineStyle','-', 'Color', 'b');
    set(H2,'LineStyle',':', 'Color', 'r');
    set(AX(1), 'XLim', [min(lambda_vec/lambda_c), max(lambda_vec/lambda_c)]);
    set(AX(2), 'XLim', [min(lambda_vec/lambda_c), max(lambda_vec/lambda_c)]);
    set(AX(1), 'YColor', 'k');
    set(AX(2), 'YColor', 'k');
    %plot(lambda_c./1e-9*ones(1,length(r)), abs(r), ':r');
    %text(lambda_c./1e-9, max(abs(r))-0.1, ['\lambda_c=' num2str(lambda_c/1e-9) ' nm']);
    %axis([min(lambda_vec/lambda_c), max(lambda_vec/lambda_c), 0, 1]);
    xlabel('\lambda/\lambda_c');
    ylabel('|r_{MC}|');
    title('(b)');
    subplot(313); hold on; box on;
    plot(phi_vec/pi, abs(r));
    %text(0.5, max(abs(r))-0.1, '\phi=\pi/2');
    axis([min(phi_vec/pi), max(phi_vec/pi), min(abs(r)), 1]);
    xlabel('\phi/\pi');
    ylabel('|r_{MC}|');
    title('(c)');
    
    figure(2);
    [maxs_r, mins_r] = peakdet(abs(r), 0.1, E_vec/Consts.e_0);

    subplot(length(con_vec),1,con_num);
    plot(E_vec/Consts.e_0, abs(r));
    axis([min(E_vec)/Consts.e_0, max(E_vec)/Consts.e_0, 0, 1]);
    for (mn=1:length(mins_r(:,1)))
        text(mins_r(mn,1), mins_r(mn,2), [num2str(mins_r(mn,1)) 'eV']);
    end
    con_text = strrep([num2str(con_vec(con_num), '%1.0e')], 'e+0', 'x10^{');
    con_text = [con_text '}'];
    ylabel(con_text);
    h_curr = gca;
    if (con_num==1)
        title_text = ['|r_{MC}|, E_c=' num2str(E_c/Consts.e_0) 'eV'];
        title(title_text);
    end
    if (con_num~=length(con_vec))
        set(h_curr, 'XTickLabel', '');
    end
    drawnow;
end
figure(2);
xlabel('E [eV]');