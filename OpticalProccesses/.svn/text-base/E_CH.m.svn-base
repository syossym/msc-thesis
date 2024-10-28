function E_CH = CalculateECH(Gq_hh,q_eps_fit,k,eps_r)

global Consts;
int_E_CH = zeros(1,length(q));
for (qq = 1:length(q))
    int_E_CH(qq) = interp1(k,Gq_hh,q(qq),'pchip').*((1./(interp1(k,eps_q,q(qq),'pchip')./q(qq)))-1);
end
D_E_CH = Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)*trapz(q,int_E_CH);
