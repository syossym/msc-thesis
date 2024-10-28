function [E_CH] = CalculateECH(Gq_hh,q_eps_fit,k,eps_r)

global Consts;

q_vec = interp(k, 10);
q_vec(1) = q_vec(2);

Gq_hh_interp = interp1(k, Gq_hh, q_vec, 'pchip');
q_eps_interp = polyval(q_eps_fit, q_vec);

C = Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0);
E_CH = C.*trapz(q_vec, Gq_hh_interp.*((1./(q_eps_interp./q_vec))-1));
