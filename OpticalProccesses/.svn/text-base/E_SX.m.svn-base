%
% This function calculates the screened-exchange shift energy
%
function D_E_SX = E_SX(Gq_ee,Gq_hh,theta,eps_q,k,f_e,f_h,eps_r)

global Consts;
C_E_SX = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
D_E_SX = zeros(1,length(k)); int_E_SX = zeros(1,length(k));
for (kk_tag = 1:length(k))
    for (kk = 1:length(k))
        k_minus_k_tag = sqrt(k(kk_tag).^2+k(kk).^2-2.*k(kk_tag).*k(kk).*cos(theta));
        form_e = interp1(k,Gq_ee,k_minus_k_tag,'pchip');
        form_e(isnan(form_e)) = Gq_ee(kk_tag);
        form_h = interp1(k,Gq_hh,k_minus_k_tag,'pchip');
        form_h(isnan(form_h)) = Gq_hh(kk_tag);
        eps = interp1(k,eps_q,k_minus_k_tag,'pchip');
        eps(isnan(eps)) = eps_q(kk_tag);
        integ = (form_e.*f_e(kk)+form_h.*f_h(kk)).*(1./eps);
%         figure(2); 
%         subplot(311); plot(theta, form_e, 'b', theta, form_h, 'r');
%         subplot(312); plot(theta, eps);
%         subplot(313); plot(theta, integ);
%         drawnow;
        int_E_SX(kk) = trapz(theta, integ);
    end
    int_E_SX = RepairNanValues(int_E_SX);
    %figure(3); plot(k, int_E_SX); drawnow;
    D_E_SX(kk_tag) = -C_E_SX.*trapz(k, k.*int_E_SX);
end

