%
% This function calculates the 2D form factor
%
function G = G_q(Xi12,k_vec,z_grid)

global z_exp;

for (kk = 1:length(k_vec))
    z_exp_q = exp(z_exp.*k_vec(kk));
    integrand_G = Xi12.*z_exp_q;
    G(kk) = trapz(z_grid, trapz(z_grid,integrand_G,1), 2);
end

%
% This function calculates the screened dielectric function
%
function eps_q = Eps_q(k_vec,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,eps_r)

global Consts;

k_len = length(k_vec);

for (kk = 1:k_len)
    int3 = zeros(size(k_vec)); int4 = int3;
    for (kk_tag = 1:k_len)
        k_tag(kk_tag,:) = sqrt(k_vec(kk_tag)^2+k_vec(kk)^2-2*k_vec(kk_tag)*k_vec(kk)*cos(theta));
        
        E_e_k_minus_q(kk_tag, :) = interp1(k_vec, E_e, k_tag(kk_tag, :));
        E_h_k_minus_q(kk_tag, :) = interp1(k_vec, E_h, k_tag(kk_tag, :));
        E_h_k_minus_q(kk_tag, isnan(E_h_k_minus_q(kk_tag, :))) = E_h(kk);
        E_e_k_minus_q(kk_tag, isnan(E_e_k_minus_q(kk_tag, :))) = E_e(kk);
        f_h_k_minus_q(kk_tag, :) = FermiDirac(E_fv, E_h_k_minus_q(kk_tag, :), T);
        f_e_k_minus_q(kk_tag, :) = FermiDirac(E_fc, E_e_k_minus_q(kk_tag, :), T);
        nom_e(kk_tag, :) = f_e_k_minus_q(kk_tag, :) - f_e(kk_tag).*ones(size(theta));
        nom_h(kk_tag, :) = f_h_k_minus_q(kk_tag, :) - f_h(kk_tag).*ones(size(theta));
        den_e(kk_tag, :) = E_e_k_minus_q(kk_tag, :) - E_e(kk_tag).*ones(size(theta));
        den_h(kk_tag, :) = E_h_k_minus_q(kk_tag, :) - E_h(kk_tag).*ones(size(theta));
        
        int1(kk_tag, :) = nom_e(kk_tag, :)./den_e(kk_tag, :);
        int2(kk_tag, :) = nom_h(kk_tag, :)./den_h(kk_tag, :);     
        int1(kk_tag, isnan(int1(kk_tag,:))) = 0; int2(kk_tag, isnan(int2(kk_tag,:))) = 0;
        
        p1 = polyfit(theta, int1(kk_tag,:), 25);
        p2 = polyfit(theta, int2(kk_tag,:), 25);
        
        int3(kk_tag) = k_vec(kk_tag).*trapz(theta, polyval(p1,theta));
        int4(kk_tag) = k_vec(kk_tag).*trapz(theta, polyval(p2,theta));
    end
    
    int3(isnan(int3))=0; int4(isnan(int4))=0;
    
    p3 = polyfit(k_vec, int3, 20);
    p4 = polyfit(k_vec, int4, 20);
    
%     figure(5);
%     subplot(211); plot(k_vec, int3, 'r.', k_vec, polyval(p3,k_vec), 'b');
%     subplot(212); plot(k_vec, int4, 'r.', k_vec, polyval(p4,k_vec), 'b');
%     drawnow;
    
    int3 =  polyval(p3,k_vec);  int4 =  polyval(p4,k_vec);
    alpha_q(kk) = (Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)).*(trapz(k_vec, int3) + trapz(k_vec, int4));
end

p_alpha = polyfit(k_vec, alpha_q, 20);
eps_q = k_vec - polyval(p_alpha,k_vec);

%
% This function calculates the screened-exchange shift energy
%
function D_E_SX = E_SX(Gq_ee,Gq_hh,theta,eps_q,z_grid,k,q,f_e,f_h,eps_r)

global Consts;
C_E_SX = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
D_E_SX = zeros(1,length(k)); int_E_SX = zeros(1,length(k));
for (kk_tag = 1:length(k))
    for (kk = 1:length(k))
        k_minus_k_tag = sqrt(k(kk_tag).^2+k(kk).^2-2.*k(kk_tag).*k(kk).*cos(theta));
        int_E_SX(kk) = trapz(theta, (interp1(k,Gq_ee,k_minus_k_tag).*f_e(kk)+...
                                     interp1(k,Gq_hh,k_minus_k_tag).*f_h(kk)).*...
                                     (1./(interp1(k,eps_q,k_minus_k_tag))));
    end
    D_E_SX(kk_tag) = -C_E_SX.*trapz(k, k.*int_E_SX);
end

%
% This function calculates the Coulomb-hole self-energy
%
function D_E_CH = E_CH(Gq_hh,eps_q,z_grid,k,q,f_e,f_h,eps_r)

global Consts;
int_E_CH = zeros(1,length(q));
for (qq = 1:length(q))
    int_E_CH(qq) = interp1(k,Gq_hh,q(qq)).*((1./(interp1(k,eps_q,q(qq))./q(qq)))-1);
end
D_E_CH = Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)*trapz(q,int_E_CH);
