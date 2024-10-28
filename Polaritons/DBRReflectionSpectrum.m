%clc; clear all; close all;

% Definitions
c = 2.99792458e8;                    % [m/s]
lambda_vec = [1300:0.1:1600].*1e-9;  % [m]
lambda_c = min(lambda_vec)+0.5*(max(lambda_vec)-min(lambda_vec))      % [m]
n_H = 3.45;
n_L = 2.98;
n_c_l = 1; %2.98;  % refractive index of the left cladding
n_c_r = 3.59;  % refractive index of the right cladding
N = 35;     % number of alternating index pairs
version = 1;

% Parameter calculation
l_L = lambda_c/n_L/4;
l_H = lambda_c/n_H/4;
n_vec = [repmat([n_H, n_L],1,N), n_H];
l_vec = [repmat([l_H, l_L],1,N), l_H];
n_profile = [n_c_l.*ones(1,30), repmat([n_H*ones(1,30), n_L*ones(1,30)],1,N), n_H*ones(1,30), n_c_r.*ones(1,30)];
dl = sum(l_vec)./length(n_profile);      % [m]
l_profile = [30*dl*(-1:1/30:0), (length(n_profile)-30)*dl*(0:1/(length(n_profile)-30):1)];

% Simulation
for (gg=1:length(lambda_vec))
    
    switch version
        case 1,
            
            M_DBR = eye(2);
            M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
            M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
            M_DBR = M_c_r*M_DBR;
            phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
            %phi_vec(gg) = 0;
            
            for (nn=length(n_vec):-1:1)
                k = (2*pi/lambda_vec(gg))*n_vec(nn);
                %phi_vec(gg) = phi_vec(gg) + (2*pi/lambda_vec(gg))*n_vec(nn)*l_vec(nn);                
                M_p = [exp(1i*k*l_vec(nn)), 0 ;...
                    0 , exp(-1i*k*l_vec(nn))];
                if (nn~=1)
                    M_i = (1/(2*n_vec(nn))).*[(n_vec(nn)+n_vec(nn-1)) , (n_vec(nn)-n_vec(nn-1)) ; ...
                        (n_vec(nn)-n_vec(nn-1)) , (n_vec(nn)+n_vec(nn-1))];
                else
                    M_i = eye(2);
                end
                
                M_DBR = M_i*M_p*M_DBR;
                
            end
            M_DBR = M_c_l*M_DBR;
            r(gg) = -M_DBR(2,1)/M_DBR(2,2);
            
        case 2,
            
            M_r = [1, 1; n_c_r, -n_c_r];
            M_l = [1, 1; n_c_l, -n_c_l];
            M_i_H = [1, 1; n_H, -n_H];
            M_i_L = [1, 1; n_L, -n_L];
            k_H = (2*pi/lambda_vec(gg))*n_H;
            k_L = (2*pi/lambda_vec(gg))*n_L;
            M_p_H = [exp(1i*k_H*l_H), 0; 0, exp(-1i*k_H*l_H)];
            M_p_L = [exp(1i*k_L*l_L), 0; 0, exp(-1i*k_L*l_L)];
            
            M_stack = M_i_H*M_p_H*inv(M_i_H)*M_i_L*M_p_L*inv(M_i_L);
            M_DBR = eye(2);
            for (nn=1:N)
                M_DBR = M_DBR*M_stack;
            end
            
            M_DBR = inv(M_l)*M_DBR*M_r;
            
            r(gg) = M_DBR(2,1)/M_DBR(1,1);
            
        case 3,
            
            M_DBR = eye(2);
            for (nn=1:2*N)
                k = (2*pi/lambda_vec(gg))*n_vec(nn);
                M_i = [cos(k*l_vec(nn)), -1i*(1/n_vec(nn))*sin(k*l_vec(nn)); ...
                    -1i*(n_vec(nn))*sin(k*l_vec(nn)), cos(k*l_vec(nn))];
                M_DBR = M_DBR*M_i;
            end
            M = M_DBR;
            r(gg) = (M(2,1)+n_c_r*M(2,2)-M(1,1)-n_c_r*M(1,2))/(M(2,1)+n_c_r*M(2,2)+M(1,1)+n_c_r*M(1,2));
    end
    
end

r_anal = 1-2*(n_c_l/n_c_r)*(n_L/n_H)^(2*N);

% Plotting
figure(1);
subplot(311);
plot(l_profile(1:length(n_profile))./1e-9, n_profile); hold on;
vec = 0:max(n_profile)/length(n_profile):max(n_profile)-max(n_profile)/length(n_profile);
plot(zeros(1,length(vec)), vec, ':g',...
     ones(1,length(vec)).*(max(l_profile)-30*dl)/1e-9, vec, ':g'); hold on;
axis([min(l_profile)./1e-9, max(l_profile)./1e-9, 0 , max(n_profile)+0.2]);
xlabel('z [nm]');
ylabel('n'); 
title('(a)');
subplot(312); 
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
ylabel('|r_{DBR}|');
title('(b)');
subplot(313); 
plot(phi_vec/pi, abs(r)); hold on; 
plot(0.5*ones(1,length(r)), abs(r), 'r');
text(0.5, max(abs(r))-0.1, '\phi=\pi/2');
axis([min(phi_vec/pi), max(phi_vec/pi), min(abs(r)), 1]);
xlabel('\phi/\pi');
ylabel('|r_{DBR}|');
title('(c)');


% % Consts
% lambda_vec = 1300:0.1:1800;       % (nm)
% lambda0 = 1550;                   % (nm)
% n0 = 3.05; %1.45;
% d0 = 0.25*(lambda0/n0);           % (nm)
% dn_vec = 3.6-3.05;  % [0.15,0.1,0.05,0.02];
% n1_vec = n0 + dn_vec;
% d1_vec = 0.25*(lambda0./n1_vec);  % (nm)
% N = 35;
%
% % Calcs
% hold on;
% for (ii = 1:length(dn_vec))
%
%     for (jj=1:length(lambda_vec))
%
%         % Consts
%         r11 = (n1_vec(ii)-n0)/(n1_vec(ii)+n0);
%         r12 = r11;
%         t11 = sqrt(1-r11^2);
%         t12 = sqrt(1-r12^2);
%         beta1_vec = (2*pi*n1_vec(ii))./lambda_vec(jj);  % (rad/nm)
%
%         r21 = (n0-n1_vec(ii))/(n1_vec(ii)+n0);
%         r22 = r21;
%         t21 = sqrt(1-r21^2);
%         t22 = sqrt(1-r22^2);
%         beta2_vec = (2*pi*n0)./lambda_vec(jj);          % (rad/nm)
%
%         r31 = (n0-n1_vec(ii))/(n1_vec(ii)+n0);
%         r32 = (n0-1)/(1+n0);
%         t31 = sqrt(1-r31^2);
%         t32 = sqrt(1-r32^2);
%         beta3_vec = (2*pi*n0)./lambda_vec(jj);          % (rad/nm)
%
%         M1_mat = [ (1/(t11*t12)).*( exp(i*beta1_vec*d1_vec(ii)) - r11*r12*exp(-i*beta1_vec*d1_vec(ii)) ), ...
%             (-1/(t11*t12)).*( r11*exp(-i*beta1_vec*d1_vec(ii)) - r12*exp(i*beta1_vec*d1_vec(ii)) ); ...
%             (-1/(t11*t12)).*( r11*exp(i*beta1_vec*d1_vec(ii)) - r12*exp(-i*beta1_vec*d1_vec(ii)) ), ...
%             (1/(t11*t12)).*( exp(-i*beta1_vec*d1_vec(ii)) - r11*r12*exp(i*beta1_vec*d1_vec(ii)) ) ];
%
%         M2_mat = [ (1/(t21*t22)).*( exp(i*beta2_vec*d0) - r21*r22*exp(-i*beta2_vec*d0) ), ...
%             (-1/(t21*t22)).*( r21*exp(-i*beta2_vec*d0) - r22*exp(i*beta2_vec*d0) ); ...
%             (-1/(t21*t22)).*( r21*exp(i*beta2_vec*d0) - r22*exp(-i*beta2_vec*d0) ), ...
%             (1/(t21*t22)).*( exp(-i*beta2_vec*d0) - r21*r22*exp(i*beta2_vec*d0) ) ];
%
%         M3_mat = [ (1/(t31*t32)).*( exp(i*beta3_vec*d0) - r31*r32*exp(-i*beta3_vec*d0) ), ...
%             (-1/(t31*t32)).*( r31*exp(-i*beta3_vec*d0) - r32*exp(i*beta3_vec*d0) ); ...
%             (-1/(t31*t32)).*( r31*exp(i*beta3_vec*d0) - r32*exp(-i*beta3_vec*d0) ), ...
%             (1/(t31*t32)).*( exp(-i*beta3_vec*d0) - r31*r32*exp(i*beta3_vec*d0) ) ];
%
%         % Stacking
%         M_mat_sec = M1_mat*M2_mat;
%         M_mat_end_sec = M1_mat*M3_mat;
%         Mtot = M_mat_end_sec;
%         for (kk=1:N-1)
%             Mtot = M_mat_sec*Mtot;
%         end
%
%         R12(ii,jj) = abs(Mtot(2,1)./Mtot(1,1)).^2;
%
%     end
%
%     plot(lambda_vec, R12(ii,:));
%     gtext([ '{\Delta}n=', num2str(dn_vec(ii))]);
%
% end
%
% xlabel( '{\lambda}, nm' );
% ylabel('Reflectivity');
% title( [' R_1_2-{\lambda} for N_s_e_c=', num2str(N)] );
% hold off;