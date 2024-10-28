function Theta = ThetaMat(Gq,theta,eps_q,k_vec,eps_r)

global Consts;

C_Theta = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
Theta = zeros(length(k_vec));
for (kk = 1:length(k_vec))
    for (kk_tag = 1:length(k_vec))
        %disp(['i=' num2str(ii) ', j=' num2str(jj)]);
        q_tag = sqrt(k_vec(kk).^2 + k_vec(kk_tag).^2 - 2.* k_vec(kk).*k_vec(kk_tag).*cos(theta));
        %q(q==0) = 1e-5;
        
        form_eh = interp1(k_vec,Gq,q_tag,'pchip');
        form_eh(isnan(form_eh)) = Gq(kk);
        eps = interp1(k_vec,eps_q,q_tag,'pchip');
        eps(isnan(eps)) = eps_q(kk);
        
        Theta(kk,kk_tag) = C_Theta.*trapz(theta, form_eh./eps);
    end
    %figure(2); plot(k_vec, Theta(kk,:), 'r');
    drawnow;
end