function [q_eps, q_eps_fit] = Screening_Test1(k_vec,theta,E_fv,E_fc,T,E_e,E_h,eps_r)

global Consts;

k_vec_interp = interp(k_vec, 10);
q_vec = interp(k_vec, 10);

E_e_interp = interp1(k_vec, E_e, k_vec_interp, 'pchip');
E_h_interp = interp1(k_vec, E_h, k_vec_interp, 'pchip');

f_e_interp = FermiDirac(E_fc, E_e_interp, T);
f_h_interp = FermiDirac(E_fv, E_h_interp, T);

E_e_fit = polyfit(k_vec_interp, E_e_interp, 6);
E_h_fit = polyfit(k_vec_interp, E_h_interp, 30);

for (qq=1:length(q_vec))
    
    for (kk=1:length(k_vec_interp))
        
        k_minus_q = sqrt(q_vec(qq)^2 + k_vec_interp(kk)^2 - 2*q_vec(qq)*k_vec_interp(kk)*cos(theta));
        
        E_e_temp = polyval(E_e_fit, k_minus_q);
        E_h_temp = polyval(E_h_fit, k_minus_q);
        f_e_temp = FermiDirac(E_fc, E_e_temp, T);
        f_h_temp = FermiDirac(E_fv, E_h_temp, T);
        
        int_e_theta(kk) = trapz(theta, (f_e_temp-f_e_interp(kk))./(E_e_temp-E_e_interp(kk)));
        int_h_theta(kk) = trapz(theta, (f_h_temp-f_h_interp(kk))./(E_h_temp-E_h_interp(kk)));
        
    end
    
    int_e_k(qq) = trapz(k_vec_interp, k_vec_interp.*int_e_theta);
    int_h_k(qq) = trapz(k_vec_interp, k_vec_interp.*int_h_theta);
    
    int_k(qq) = int_e_k(qq) + int_h_k(qq);
    
end

q_eps = q_vec - (Consts.e_0^2/(4*pi^2*Consts.eps_0*R.Materials{2}.eps_r))*int_k;

q_eps_fit = polyfit(q_vec, q_eps, 20);

