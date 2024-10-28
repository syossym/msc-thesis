function [Theta, E_SX] = CalculateThetaESX(Gq_ee,Gq_hh,Gq_eh,theta,q_eps_fit,k,f_e,f_h,eps_r)

global Consts;

k_vec_interp = interp(k, 1);
k_tag_vec_interp = interp(k, 1);

Gq_ee_interp = interp1(k, Gq_ee, k_vec_interp, 'pchip');
Gq_hh_interp = interp1(k, Gq_hh, k_vec_interp, 'pchip');
Gq_eh_interp = interp1(k, Gq_eh, k_vec_interp, 'pchip');

Gq_ee_fit = polyfit(k_vec_interp, Gq_ee_interp, 20);
Gq_hh_fit = polyfit(k_vec_interp, Gq_hh_interp, 20);
Gq_eh_fit = polyfit(k_vec_interp, Gq_eh_interp, 20);

f_e_interp = interp1(k, f_e, k_tag_vec_interp, 'pchip');
f_h_interp = interp1(k, f_h, k_tag_vec_interp, 'pchip');

C = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
E_SX = zeros(1,length(k)); 
int_E_SX = zeros(1,length(k_vec_interp));
Theta = zeros(length(k));
q_eps_full = polyval(q_eps_fit, k_vec_interp);
eps_full = q_eps_full./k_vec_interp;
eps_full(1) = eps_full(2);

for (kk = 1:length(k_vec_interp))

    for (kk_tag = 1:length(k_tag_vec_interp))
        q = sqrt(k_vec_interp(kk).^2 + k_tag_vec_interp(kk_tag).^2 - 2.*k_vec_interp(kk).*k_tag_vec_interp(kk_tag).*cos(theta));
       
        %q(q > max(k_vec_interp)) = max(k_vec_interp);
        %q(q < min(k_vec_interp)) = min(k_vec_interp);
        
        % Theta
        form_eh = polyval(Gq_eh_fit, q);
        %q_eps = polyval(q_eps_fit, q);
        q_eps = interp1(k_vec_interp, q_eps_full, q, 'pchip');
        eps = interp1(k_vec_interp, eps_full, q, 'pchip');
        Theta_1(kk,kk_tag) = C.*k_tag_vec_interp(kk_tag).*trapz(theta, form_eh./q_eps);
        Test(kk,kk_tag) = C.*k_tag_vec_interp(kk_tag).*trapz(theta, 1./q);
        
%         figure(1);
%         subplot(211); plot(k, Gq_eh, 'b', q, form_eh, 'r.'); 
%         subplot(212); plot(k, polyval(q_eps_fit, k), 'b', q, q_eps, 'r.'); drawnow;
        
        % E_SX
        form_ee = polyval(Gq_ee_fit, q);
        form_hh = polyval(Gq_hh_fit, q);
        int_E_SX(kk_tag) = trapz(theta, (form_ee.*f_e_interp(kk_tag)+form_hh.*f_h_interp(kk_tag))./q_eps); 
    end
    E_SX_1(kk) = -C.*trapz(k_tag_vec_interp, k_tag_vec_interp.*int_E_SX);
end

for (kk = 1:length(k_vec_interp))
   Theta_2(kk,:) = interp1(k_tag_vec_interp, Theta_1(kk,:), k, 'pchip');
end

for (kk = 1:length(Theta_2(1,:)))
   Theta(:,kk) = interp1(k_vec_interp.', Theta_2(:,kk), k', 'pchip');
end

E_SX = interp1(k_vec_interp, E_SX_1, k, 'pchip');