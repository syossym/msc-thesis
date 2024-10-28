function [Theta] = CalculateTheta(Gq_eh,q_eps_fit,k,eps_r)

global Consts;

k_vec_interp = interp(k, 1);

Gq_eh_interp = interp1(k, Gq_eh, k_vec_interp, 'pchip');
Gq_eh_fit = polyfit(k_vec_interp, Gq_eh_interp, 20);

f_e_interp = interp1(k, f_e, k_tag_vec_interp, 'pchip');
f_h_interp = interp1(k, f_h, k_tag_vec_interp, 'pchip');

C = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
Theta = zeros(length(k));
q_eps_full = polyval(q_eps_fit, k_vec_interp);
eps_full = q_eps_full./k_vec_interp;
eps_full(1) = eps_full(2);

for (kk = 1:length(k_vec_interp))

    for (kk_tag = 1:length(k_tag_vec_interp))
        q = sqrt(k_vec_interp(kk).^2 + k_tag_vec_interp(kk_tag).^2 - 2.*k_vec_interp(kk).*k_tag_vec_interp(kk_tag).*cos(theta));
       
        %q(q > max(k_vec_interp)) = max(k_vec_interp);
        %q(q < min(k_vec_interp)) = min(k_vec_interp);

        form_eh = polyval(Gq_eh_fit, q);
        %q_eps = polyval(q_eps_fit, q);
        q_eps = interp1(k_vec_interp, q_eps_full, q, 'pchip');
        eps = interp1(k_vec_interp, eps_full, q, 'pchip');
        Theta_1(kk,kk_tag) = C.*trapz(theta, 1./q_eps);
        Test(kk,kk_tag) = C.*trapz(theta, 1./q);
        
%         figure(1);
%         subplot(211); plot(k, Gq_eh, 'b', q, form_eh, 'r.'); 
%         subplot(212); plot(k, polyval(q_eps_fit, k), 'b', q, q_eps, 'r.'); drawnow;
        
    end
end

for (kk = 1:length(k_vec_interp))
   Theta_2(kk,:) = interp1(k_tag_vec_interp, Theta_1(kk,:), k, 'pchip');
end

for (kk = 1:length(Theta_2(1,:)))
   Theta(:,kk) = interp1(k_vec_interp.', Theta_2(:,kk), k', 'pchip');
end
