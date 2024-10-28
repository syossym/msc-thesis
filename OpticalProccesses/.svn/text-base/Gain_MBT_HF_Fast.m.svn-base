function [G_TE, G_TM, R_sp, delta_n, deltaE_ch] = Gain_MBT_HF_Fast(E_c, E_v, Efc, Efv, cb_num, vb_num, k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, q_min, dq, q_point, omega, omega_point, Structure, G_mat, q_times_epsilon_q, ukE, ukM, T, gamma)

global Consts;

UT = Consts.k_B*T;
uk0_TE = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
uk0_TM = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
Eshift = 0;
Epsr_static = Structure{2}.eps_r;%*Consts.eps_0;
n_opt = Structure{2}.n;
L = Structure{2}.Width*1e-10;

k_normal = zeros(1, k_point);
dk_normal = zeros(1, k_point);
kk = zeros(1, k_point);
dkk = zeros(1, k_point);
omega_k = zeros(1, k_point);
A0 = zeros(1, k_point);
B0 = zeros(1, k_point);
ukE0 = zeros(1, k_point);
ukM0 = zeros(1, k_point);
deltaE_ch = zeros(1, vb_num);
deltaE_sx_e0 = zeros(cb_num, k_point);
deltaE_sx_h0 = zeros(vb_num, k_point);
deltaE_sx_e = zeros(1, k_point);
deltaE_sx_h = zeros(1, k_point);
theta = zeros(k_point, k_prime_point);
theta0 = zeros(k_point, k_prime_point);
theta_help = zeros(k_point, k_prime_point);

v2 = zeros(1,k_point);
v3 = zeros(1,k_point);

for (I=1:k_point)
    k_normal(I) = k_min+(I-1)*dk;
    if ((I==1)||(I==k_point))
        dk_normal(I) = dk/2;
    else
        dk_normal(I) = dk;
    end
end

for (J1=1:cb_num)
    
    deltaE_sx_e0(J1,:) = FindDeltaEsx(k_normal, dk_normal, k_point, k_normal, dk_normal, k_point, E_c(J1,:), Efc, G_mat{J1,J1}, q_times_epsilon_q, dq, Epsr_static, T);
    
    for (J2=1:vb_num)
        [J1,J2]
        
        if (J1==1)
            deltaE_ch(J2) = FindDeltaEch(q_min, dq, q_point, q_times_epsilon_q, Epsr_static, G_mat{J2+cb_num, J2+cb_num});
            deltaE_sx_h0(J2,:) = FindDeltaEsx(k_normal, dk_normal, k_point, k_normal, dk_normal, k_point, E_v(J2,:), Efv, G_mat{J2+cb_num,J2+cb_num}, q_times_epsilon_q, dq, Epsr_static, T);
        end
        
        Eshift = Eshift + (deltaE_sx_e0(J1,1) + deltaE_sx_h0(J2,1) + deltaE_ch(J2))/(cb_num*vb_num);
        theta0 = FindTheta(k_normal, k_point, k_normal, k_point, G_mat{J1,cb_num+J2}, q_times_epsilon_q, dq, Epsr_static);
        
        for (nw=1:k_point)
            omega_k(nw) = (E_c(J1,nw) + E_v(J2,nw) + deltaE_ch(J2) + deltaE_sx_e0(J1,nw) + deltaE_sx_h0(J2,nw))/Consts.hbar;
        end
        
        omega_k_min = omega_k(1);
        omega_k_max = omega_k(k_point);
        
        for (R=1:omega_point)
            
            C0 = 1/(omega(R)*Consts.eps_0*n_opt*Consts.c*Consts.hbar*pi*L*2)/100;      % 100 because m-1->cm-1
            
            if ((omega(R)<=omega_k_min)||(omega(R)>=omega_k_max))
                for (nk=1:k_point)
                    kk(nk) = k_normal(nk);
                    dkk(nk) = dk_normal(nk);
                    A0(nk) = E_c(J1,nk);
                    B0(nk) = E_v(J2,nk);
                    ukE0(nk) = squeeze(ukE(J1,J2,nk));
                    ukM0(nk) = squeeze(ukM(J1,J2,nk));
                    deltaE_sx_e(nk) = deltaE_sx_e0(J1,nk);
                    deltaE_sx_h(nk) = deltaE_sx_h0(J2,nk);
                    for (nk_prime=1:k_point)
                        theta(nk,nk_prime) = theta0(nk,nk_prime);
                    end
                end
            else
                
                kk = FindK(k_point, k_normal(1), k_normal(k_point), omega(R), omega_k, k_normal);
                
                for (nk=1:k_point)
                    if (nk<k_point)
                        v1 = kk(nk)/dk;
                        v2(nk) = floor(v1);
                        v3(nk) = v1-v2(nk);
                    else
                        v2(k_point) = k_point-2;
                        v3(k_point) = 1;
                    end
                    
                    if (nk==1)
                        dkk(1) = (kk(2)-kk(1))/2;
                    else
                        if (nk==k_point)
                            dkk(k_point) = (kk(k_point)-kk(k_point-1))/2;
                        else
                            dkk(nk) = (kk(nk+1)-kk(nk-1))/2;
                        end
                    end
                    
                    if (v2(nk)+2 <= k_point)
                        A0(nk) = (1-v3(nk))*E_c(J1,v2(nk)+1) + v3(nk)*E_c(J1,v2(nk)+2);
                        B0(nk) = (1-v3(nk))*E_v(J2,v2(nk)+1) + v3(nk)*E_v(J2,v2(nk)+2);
                        ukE0(nk) = (1-v3(nk))*squeeze(ukE(J1,J2,v2(nk)+1)) + v3(nk)*squeeze(ukE(J1,J2,v2(nk)+2));
                        ukM0(nk) = (1-v3(nk))*squeeze(ukM(J1,J2,v2(nk)+1)) + v3(nk)*squeeze(ukM(J1,J2,v2(nk)+2));
                        deltaE_sx_e(nk) = (1-v3(nk))*deltaE_sx_e0(J1,v2(nk)+1) + v3(nk)*deltaE_sx_e0(J1,v2(nk)+2);
                        deltaE_sx_h(nk) = (1-v3(nk))*deltaE_sx_h0(J2,v2(nk)+1) + v3(nk)*deltaE_sx_h0(J1,v2(nk)+2);
                    else
                        A0(nk) = A0(nk-1);
                        B0(nk) = B0(nk-1);
                        ukE0(nk) = ukE0(nk-1);
                        ukM0(nk) = ukM0(nk-1);
                        deltaE_sx_e(nk) = deltaE_sx_e(nk-1);
                        deltaE_sx_h(nk) = deltaE_sx_h(nk-1);
                    end
                end
                
                for (Nk=1:k_point)
                    for (Nk_prime=1:k_point)
                        if (v2(Nk_prime)+2 <= k_point)
                            theta(Nk,Nk_prime) = (1-v3(Nk_prime))*theta0(Nk,v2(Nk_prime)+1) + v3(Nk_prime)*theta0(Nk,v2(Nk_prime)+2);
                        else
                            theta(Nk,Nk_prime) = theta(Nk,Nk_prime-1);
                        end
                    end
                end
                
                for (Nk=1:k_point)
                    for (Nk_prime=1:k_point)
                        if (v2(Nk_prime)+2 <= k_point)
                            theta_help(Nk_prime) = (1-v3(Nk_prime))*theta(v2(Nk_prime)+1,Nk) + v3(Nk_prime)*theta(v2(Nk_prime)+2,Nk);
                        else
                            theta_help(Nk_prime) = theta_help(Nk_prime-1);
                        end
                    end
                    for (Nk_prime=1:k_point)
                        theta(Nk_prime,Nk) = theta_help(Nk_prime);
                    end
                end
            end
            
            [QkE, QkM, Qksp] = Enhancement(theta, kk, dkk, k_point, kk, dkk, k_point, A0, Efc, B0, Efv, T, omega(R), ukE0, ukM0, deltaE_ch(J2), deltaE_sx_e, deltaE_sx_h, uk0_TE, uk0_TM, gamma, 'Lorentzian_enhanced');
            
            if ((J1==1)&&(J2==1))
                G_TE(R) = 0;
                G_TM(R) = 0;
                R_sp(R) = 0;
            end
            
            for (I=1:k_point)
                G_TE(R) = G_TE(R) + 2*C0*real(QkE(I))*2*dkk(I);
                G_TM(R) = G_TM(R) + 2*C0*real(QkM(I))*2*dkk(I);
                R_sp(R) = R_sp(R) + 2*C0*real(Qksp(I))*2*dkk(I);
                delta_n(R) = 0;    % need more k-points to calculate refractive index change
            end
        end
    end
end

%%% -------------------------------------------------------------------------------------------------------------------------------%%%

function deltaE_sx = FindDeltaEsx(k, dk, k_point, k_prime, dk_prime, k_prime_point, A_B, Ef, G, q_times_epsilon_q, dq, Epsr, T)

global Consts;

UT = Consts.k_B*T; %/Consts.e_0;
Co = (Consts.e_0/pi)^2/8/Consts.eps_0/Epsr;

theta = zeros(k_point, k_prime_point);
deltaE_sx = zeros(1,k_point);
Fermi = zeros(1,k_prime_point);

phi_min = 0;
phi_max = pi;
phi_point = 51;
dphi = (phi_max-phi_min)/(phi_point-1);

for (I=1:k_point)
    deltaE_sx(I) = 0;
    for (J=1:k_prime_point)
        if (J<I)
            theta(I,J) = theta(J,I);
        else
            theta(I,J) = 0;
            if (I==1)
                Fermi(J) = 1/(1 + exp((A_B(J)-Ef)/UT));
            end
            for (K=1:phi_point)
                phi = phi_min + (K-1)*dphi;
                q_prime = sqrt(abs(k(I)^2+k_prime(J)^2-2*k(I)*k_prime(J)*cos(phi)));
                v1 = q_prime/dq;
                v2 = floor(v1);
                v3 = v1-v2;
                if (v2+2<length(k))
                    F = ((1-v3)*G(v2+1) + v3*G(v2+2))/((1-v3)*q_times_epsilon_q(v2+1) + v3*q_times_epsilon_q(v2+2));
                else
                    F = 0;
                end
                if ((K==1)||(K==phi_point))
                    theta(I,J) = theta(I,J) + Co*F*(2*dphi)/2;
                else
                    theta(I,J) = theta(I,J) + Co*F*(2*dphi);
                end               
            end
        end
        if ((J==1)||(J==k_prime_point))
            deltaE_sx(I) = deltaE_sx(I) - k_prime(J)*Fermi(J)*theta(I,J)*(dk_prime(J)/2); %/Consts.e_0;
        else
            deltaE_sx(I) = deltaE_sx(I) - k_prime(J)*Fermi(J)*theta(I,J)*dk_prime(J); %/Consts.e_0;
        end
    end
    
end

function deltaE_ch = FindDeltaEch(q_min, dq, point, q_times_epsilon_q, Epsr, G)

global Consts;

deltaE_ch = 0;
Co = Consts.e_0^2/4/pi/Consts.eps_0/Epsr;
dq = dq/50;
point = 50*(point-1);

for (I=1:point)
    q = q_min  + (I-1)*dq;
    v1 = q/(50*dq);
    v2 = floor(v1);
    v3 = v1-v2;
    if ((I==1)||(I==point))
        deltaE_ch = deltaE_ch + Co*((1-v3)*G(v2+1)+v3*G(v2+2))*(q/((1-v3)*q_times_epsilon_q(v2+1)+v3*q_times_epsilon_q(v2+2))-1)*(dq/2);
    else
        deltaE_ch = deltaE_ch + Co*((1-v3)*G(v2+1)+v3*G(v2+2))*(q/((1-v3)*q_times_epsilon_q(v2+1)+v3*q_times_epsilon_q(v2+2))-1)*dq;
    end
end

%deltaE_ch = deltaE_ch/Consts.e_0;

function theta = FindTheta(k, k_point, k_prime, k_prime_point, G, q_times_epsilon_q, dq, Epsr_static)

global Consts;

theta = zeros(k, k_prime);

phi_min = 0;
phi_max = pi;
phi_point = 51;
dphi=(phi_max-phi_min)/(phi_point-1);
Co = (Consts.e_0/pi)^2/8/Consts.eps_0/Epsr_static;

for (I=1:k_point)
    for (J=1:k_prime_point)
        if (J<I)
            theta(I,J) = theta(J,I);
        else
            theta(I,J) = 0;
            for (K=1:phi_point)
                phi = phi_min + (K-1)*dphi;
                q_prime2 = k(I)^2 + k_prime(J)^2 - 2*k(I)*k_prime(J)*cos(phi);
                if (q_prime2<0)
                    q_prime = 0;
                else
                    q_prime = sqrt(q_prime2);
                end
                v1 = q_prime/dq;
                v2 = floor(v1);
                v3 = v1-v2;
                if (v2+2 < length(k))
                    F = ((1-v3)*G(v2+1) + v3*G(v2+2))/((1-v3)*q_times_epsilon_q(v2+1) + v3*q_times_epsilon_q(v2+2));
                else
                    F = 0;
                end
                if ((K==1)||(K==phi_point))
                    theta(I,J) = theta(I,J) + Co*F*(2*(dphi/2));       % phi=0->2pi => factor 2 */
                else
                    theta(I,J) = theta(I,J) + Co*F*(2*dphi);           % phi=0->2pi => factor 2 */
                end
            end
        end
    end
end

function kt = FindK(k_point, k_min, k_max, frequency, frequency_k, k_norm)

a_min = 1;
a_max = 25;
da = (a_max-a_min);
I = 1;

while (frequency_k(I)<frequency)
    I = I+1;
end

k2 = k_norm(I);
k1 = k_norm(I-1);
w2 = frequency_k(I);
w1 = frequency_k(I-1);

k0 = k1 + (k2-k1)*sqrt((frequency-w1)/(w2-w1));

if (k0<=(k_max-k_min)/2)
    Nleft = ceil((k0-k_min)/((k_max-k_min)/(k_point-1)))+2;
    Nright = k_point-Nleft;
else
    Nright = ceil((k_max-k0)/((k_max-k_min)/(k_point-1)))+2;
    Nleft = k_point-Nright;
end

for (nl=1:Nleft)
    kt(nl) = k_min + (k0-k_min)*log(a_min+nl*da/(Nleft-1))/log(a_max);
end

for(nr=Nleft:k_point)
    kt(nr) = k_max - (k_max-k0)*log(a_max-(nr-Nleft)*da/(Nright-1))/log(a_max);
end

function [QkE, QkM, Qksp] = Enhancement(theta, k, dk, k_point, k_prime, dk_prime, k_prime_point, E_c, Efc, E_v, Efv, T, omega, ukE, ukM, deltaE_ch, deltaE_sx_e, deltaE_sx_h, uk0_TE, uk0_TM, gamma, broad_type)

global Consts;

UT = Consts.k_B*T; %/Consts.e_0;
omega_k_prime = zeros(1,k_prime_point);
w_k_prime = zeros(1,k_prime_point);
QkE = zeros(1,k_point);
QkM = zeros(1,k_point);
Qksp = zeros(1,k_point);

RE = zeros(1,k_prime_point);
RM = zeros(1,k_prime_point);
Rsp = zeros(1,k_prime_point);

omega_min = (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1))/Consts.hbar;
omega_max = (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1) + 0.1*Consts.e_0)/Consts.hbar;

for (I=1:k_point)
    
    Re_M_TE = 0;
    Im_M_TE = 0;
    Re_M_TM = 0;
    Im_M_TM = 0;
    
    for (J=1:k_prime_point)
        
        if (I==1)
            omega_k_prime(J) = (E_c(J) + E_v(J) + deltaE_ch + deltaE_sx_e(J) + deltaE_sx_h(J))/Consts.hbar;
            w_k_prime(J) = 1/(1+exp((E_c(J)-Efc)/UT)) + 1/(1+exp((E_v(J)-Efv)/UT))-1;
            
            if (strcmp(broad_type,'Lorentzian_enhanced'))
                if (omega<=omega_min)
                    RE(J) = (1-exp((Consts.hbar*omega - (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1)))/UT))*...
                            (uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RM(J) = (1-exp((Consts.hbar*omega - (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1)))/UT))*...
                            (uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                elseif ((omega>omega_min) && (omega<=omega_max))
                    RE1 = (1-exp((Consts.hbar*omega - (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1)))/UT))*...
                          (uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RM1 = (1-exp((Consts.hbar*omega - (Efc + Efv + deltaE_ch + deltaE_sx_e(1) + deltaE_sx_h(1)))/UT))*...
                          (uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RE2 = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                    RM2 = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
                    RE(J) = RE1 + (RE2-RE1)*(omega-omega_min)/(omega_max-omega_min);
                    RM(J) = RM1 + (RM2-RM1)*(omega-omega_min)/(omega_max-omega_min);
                else
                    RE(J) = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                    RM(J) = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
                end
            else
                RE(J) = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                RM(J) = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
            end
            Rsp(J) = 2*(uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT))/3 + ...
                       (uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_c(J)-Efv)/UT))/3;
        end
        
        M = (theta(I,J)*k_prime(J)*w_k_prime(J)*dk_prime(J)/Consts.hbar)/((omega_k_prime(J)-omega).^2 + gamma^2);
        Re_M_TE = Re_M_TE + M*ukE(J)/ukE(I)*(omega_k_prime(J)-omega);
        Im_M_TE = Im_M_TE + M*ukE(J)/ukE(I)*gamma;
        Re_M_TM = Re_M_TM + M*ukM(J)/ukM(I)*(omega_k_prime(J)-omega);
        Im_M_TM = Im_M_TM + M*ukM(J)/ukM(I)*gamma;
    end
    
    switch (broad_type)
        case 'Lorentzian',
            lineshape = gamma^2/((omega_k_prime(I)-omega).^2+gamma^2);
        case 'Lorentzian_enhanced',
            lineshape = gamma^2/((omega_k_prime(I)-omega).^2+gamma^2);
        case 'Sech',
            lineshape = 1/cosh((omega_k_prime(I)-omega)/gamma);
        case 'Cos_Hyp',
            lineshape = 1/cosh((omega_k_prime(I)-omega)/gamma);
        otherwise,
            disp('Error');
    end
    
    QkE(I) = (RE(I)*lineshape./gamma^2)./(((1+Re_M_TE)+1i*Im_M_TE)).*(gamma - 1i*(omega_k_prime(I)-omega));
    QkM(I) = (RM(I)*lineshape./gamma^2)./(((1+Re_M_TM)+1i*Im_M_TM)).*(gamma - 1i*(omega_k_prime(I)-omega));
    Qksp(I) = (Rsp(I)*lineshape./gamma^2)./((1+2*Re_M_TE/3+Re_M_TM/3) + 1i*(2*Im_M_TE/3+Im_M_TM/3)).*(gamma-1i*(omega_k_prime(I)-omega));
    
end
