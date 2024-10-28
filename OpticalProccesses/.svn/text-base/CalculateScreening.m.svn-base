function [q_eps, q_eps_fit] = CalculateScreening(k_vec,theta,E_fv,E_fc,T,E_e,E_h,eps_r,Gq_ee,Gq_hh)

global Consts;

k_vec_interp = interp(k_vec, 1);
q_vec = interp(k_vec, 1);

for (ee=1:length(E_e))
    E_e_interp{ee} = interp1(k_vec, E_e{ee}, k_vec_interp, 'pchip');
    f_e_interp{ee} = FermiDirac(E_fc, E_e_interp{ee}, T);
    E_e_fit{ee} = polyfit(k_vec_interp, E_e_interp{ee}, 6);
    %Gq_ee_interp{ee,ee} = interp1(k_vec, Gq_ee{1,ee}, q_vec, 'pchip');
end
for (hh=1:length(E_h))
    E_h_interp{hh} = interp1(k_vec, E_h{hh}, k_vec_interp, 'pchip');
    f_h_interp{hh} = FermiDirac(E_fv, E_h_interp{hh}, T);
    E_h_fit{hh} = polyfit(k_vec_interp, E_h_interp{hh}, 30);
    %Gq_hh_interp{hh,hh} = interp1(k_vec, Gq_hh{1,hh}, q_vec, 'pchip');
end

int_k = zeros(1,length(q_vec));
for (qq=1:length(q_vec))
    
    for (kk=1:length(k_vec_interp))
        
        k_minus_q = sqrt(q_vec(qq)^2 + k_vec_interp(kk)^2 - 2*q_vec(qq)*k_vec_interp(kk)*cos(theta));
        
        k_minus_q(k_minus_q > max(k_vec)) = max(k_vec);
        k_minus_q(k_minus_q < min(k_vec)) = min(k_vec);
        
        for (ee = 1:length(E_e))
            E_e_temp = polyval(E_e_fit{ee}, k_minus_q);
            f_e_temp = FermiDirac(E_fc, E_e_temp, T);
            int_e_theta{ee}(kk) = trapz(theta, (f_e_temp-f_e_interp{ee}(kk))./(E_e_temp-E_e_interp{ee}(kk)));
        end
        for (hh = 1:length(E_h))
            E_h_temp = polyval(E_h_fit{hh}, k_minus_q);
            f_h_temp = FermiDirac(E_fv, E_h_temp, T);
            int_h_theta{hh}(kk) = trapz(theta, (f_h_temp-f_h_interp{hh}(kk))./(E_h_temp-E_h_interp{hh}(kk)));
        end
    end
    
    for (ee = 1:length(E_e))
        int_e_k{ee}(qq) = trapz(k_vec_interp, k_vec_interp.*int_e_theta{ee});
        int_k(qq) = int_k(qq) + int_e_k{ee}(qq);
    end
    for (hh = 1:length(E_h))
        int_h_k{hh}(qq) = trapz(k_vec_interp, k_vec_interp.*int_h_theta{hh});
        int_k(qq) = int_k(qq) + int_h_k{hh}(qq);
    end
   
end

q_eps = q_vec - (Consts.e_0^2/(4*pi^2*Consts.eps_0*eps_r)).*smooth(int_k, 0.1, 'rloess').';
q_eps_fit = polyfit(q_vec, q_eps, 20);
q_eps = interp1(q_vec, q_eps, k_vec, 'pchip');

