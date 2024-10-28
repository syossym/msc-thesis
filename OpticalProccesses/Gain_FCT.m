function [G_TE, G_TM, R_sp, delta_n] = Gain_FCT(E_c, E_v, Efc, Efv, cb_num, vb_num, k_min, dk, k_point, omega, omega_point, Structure, ukE, ukM, T, gamma, broad_type)

global Consts;

UT = Consts.k_B*T;
uk0_TE = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
uk0_TM = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
n_opt = Structure{2}.n;
L = Structure{2}.Width*1e-10;

omega_k = zeros(1, k_point);
G_TE = zeros(1, omega_point);
G_TM = zeros(1, omega_point);
R_sp = zeros(1, omega_point);
delta_n = zeros(1, omega_point);

for (J1=1:cb_num)
    for (J2=1:vb_num)
        for (R=1:omega_point)
            C0 = 1/(omega(R)*Consts.eps_0*n_opt*Consts.c*Consts.hbar*pi*L*2)/100;    % 100 because m-1->cm-1
            Cn = 1/(omega(R)*Consts.eps_0*n_opt*Consts.hbar*pi*L*2)/100;
            if ((J1==0)&&(J2==0))
                G_TE(R) = 0;
                G_TM(R) = 0;
                delta_n(R) = 0;
                R_sp(R) = 0;
            end
            for (I=1:k_point)
                if (R==1)
                    k = k_min+(I-1)*dk;
                    omega_k(I) = (E_c(J1,I) + E_v(J2,I))/Consts.hbar;
                    if (strcmp(broad_type, 'Lorentzian_enhanced'))
                        %RE(I) = (uk0_TE*squeeze(ukE(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT)).*(1/(1+exp((E_v(J2,I)-Efv)/UT))));
                        %RM(I) = (uk0_TM*squeeze(ukM(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT)).*(1/(1+exp((E_v(J2,I)-Efv)/UT))));
                        RE(I) = (uk0_TE*ukE{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT)).*(1/(1+exp((E_v(J2,I)-Efv)/UT))));
                        RM(I) = (uk0_TM*ukM{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT)).*(1/(1+exp((E_v(J2,I)-Efv)/UT))));
                    else
                        %RE(I) = (uk0_TE*squeeze(ukE(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))+1/(1+exp((E_v(J2,I)-Efv)/UT))-1);
                        %RM(I) = (uk0_TM*squeeze(ukM(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))+1/(1+exp((E_v(J2,I)-Efv)/UT))-1);
                        RE(I) = (uk0_TE*ukE{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))+1/(1+exp((E_v(J2,I)-Efv)/UT))-1);
                        RM(I) = (uk0_TM*ukM{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))+1/(1+exp((E_v(J2,I)-Efv)/UT))-1);                 
                    end
                    %Rsp(I) = 2*(uk0_TE*squeeze(ukE(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))*1/(1+exp((E_v(J2,I)-Efv)/UT)))/3 + ...
                    %           (uk0_TM*squeeze(ukM(J1,J2,I))).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))*1/(1+exp((E_v(J2,I)-Efv)/UT)))/3;
                    Rsp(I) = 2*(uk0_TE*ukE{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))*1/(1+exp((E_v(J2,I)-Efv)/UT)))/3 + ...
                               (uk0_TM*ukM{J1,J2}(I)).^2.*k.*(1/(1+exp((E_c(J1,I)-Efc)/UT))*1/(1+exp((E_v(J2,I)-Efv)/UT)))/3;
                
                end
                
                switch (broad_type)
                    case 'Lorentzian',
                        lineshape = (gamma^2)/((omega_k(I)-omega(R)).^2 + gamma^2);
                        lineshape_sp = lineshape;
                        lineshape_dn = lineshape;
                    case 'Lorentzian_enhanced'
                        lineshape = (1-exp((Consts.hbar*omega(R) - (Efc+Efv))/UT))*(gamma^2)/((omega_k(I)-omega(R)).^2+gamma^2);
                        lineshape_sp = gamma^2/((omega_k(I)-omega(R)).^2+gamma^2);
                        lineshape_dn = lineshape_sp;
                    case 'Sech',
                        lineshape = 1/cosh((omega_k(I)-omega(R))/gamma);
                        lineshape_sp = lineshape;
                        lineshape_dn = gamma^2/((omega_k(I)-omega(R)).^2+gamma^2);
                    case 'Cos_Hyp',
                        lineshape = 1/cosh((omega_k(I)-omega(R))/gamma);
                        lineshape_sp = lineshape;
                        lineshape_dn = gamma^2/((omega_k(I)-omega(R)).^2+gamma^2);
                end
                
                G_TE(R) = G_TE(R) + 2*C0*real(RE(I).*(gamma-1i*(omega_k(I)-omega(R))))*lineshape/gamma^2*2*dk;
                G_TM(R) = G_TM(R) + 2*C0*real(RM(I).*(gamma-1i*(omega_k(I)-omega(R))))*lineshape/gamma^2*2*dk;
                R_sp(R) = R_sp(R) + 2*C0*real(Rsp(I).*(gamma-1i*(omega_k(I)-omega(R))))*lineshape_sp/gamma^2*2*dk;
                delta_n(R) = delta_n(R) - Cn/omega(R)*RE(I)*(-(omega_k(I)-omega(R)))*lineshape_dn/gamma^2*2*dk;
            end
        end
    end
end

