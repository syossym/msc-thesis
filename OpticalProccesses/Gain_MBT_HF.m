function [G_TE, G_TM, R_sp, delta_n, deltaE_ch, deltaE_sx, Q_mat, delta_renorm] = Gain_MBT_HF(E_c, E_v, Efc, Efv, cb_num, vb_num, k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, q_min, dq, q_point, omega, omega_point, Structure, G_mat, q_times_epsilon_q, ukE, ukM, T, gamma, Theta)

global Consts;

UT = Consts.k_B*T;
uk0_TE = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
uk0_TM = sqrt(Consts.m_0*Structure{2}.E_p*Consts.e_0/6)*Consts.e_0/Consts.m_0;
Epsr_static = Structure{2}.eps_r;
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

G_TE = zeros(1, omega_point);
G_TM = zeros(1, omega_point);
R_sp = zeros(1, omega_point);
delta_n = zeros(1, omega_point);

% q_times_epsilon_q = [k_min:dk: dk*(k_point-1)];
% for (ii=1:length(E_c))
%     for (jj=1:length(E_v))
%         temp = CalculateScreening(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, E_c(ii,:), E_v(jj,:), Efc, Efv, Epsr_static, T);
%         q_times_epsilon_q = q_times_epsilon_q + temp;
%     end
% end
% 
% q_times_epsilon_q(2) = q_times_epsilon_q(1);

for (J1=1:cb_num)
    
    deltaE_sx_e(J1,:) = FindDeltaEsx(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, E_c(J1,:), Efc, G_mat{J1,J1}, q_times_epsilon_q, dq, Epsr_static, T);
    
    for (J2=1:vb_num)
        %[J1, J2]
        
        if (J1==1)
            deltaE_ch(J2) = FindDeltaEch(q_min, dq, q_point, q_times_epsilon_q, Epsr_static, G_mat{J2+cb_num,J2+cb_num});
            deltaE_sx_h(J2,:) = FindDeltaEsx(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, E_v(J2,:), Efv, G_mat{J2+cb_num,J2+cb_num}, q_times_epsilon_q, dq, Epsr_static, T);
        end
        %theta = FindTheta(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, G_mat{J1,cb_num+J2}, q_times_epsilon_q, dq, Epsr_static);
        theta = Theta{J1,J2};

        deltaE_sx{J1,J2} = deltaE_sx_e(J1,1) + deltaE_sx_h(J2,1);
        
        delta_renorm{J1,J2} = (deltaE_ch(J2).*ones(size(deltaE_sx{J1,J2})) + deltaE_sx{J1,J2});
        
        for (R=1:omega_point)
            C0 = 1/(omega(R)*Consts.eps_0*n_opt*Consts.c*Consts.hbar*pi*L*2)/100;                 % 100 because m-1->cm-1
%             [QkE, QkM, Qksp, Q] = Enhancement(theta, k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, ...
%                 E_c(J1,:), Efc, E_v(J2,:), Efv, T, omega(R), squeeze(ukE(J1,J2,:)), squeeze(ukM(J1,J2,:)), ...
%                 deltaE_ch(J2), deltaE_sx_e(J1,:), deltaE_sx_h(J2,:), uk0_TE, uk0_TM, gamma, 'Lorentzian');
            [QkE, QkM, Qksp, Q] = Enhancement(theta, k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, ...
                                              E_c(J1,:), Efc, E_v(J2,:), Efv, T, omega(R), ukE{J1,J2}, ukM{J1,J2}, ...
                                              deltaE_ch(J2), deltaE_sx_e(J1,:), deltaE_sx_h(J2,:), uk0_TE, uk0_TM, gamma, 'Lorentzian');
            Q_mat(R) = Q; 
            
            if ((J1==1)&&(J2==1))
                G_TE(R) = 0;
                G_TM(R) = 0;
                R_sp(R) = 0;
            end
            for (I=1:k_point)
                G_TE(R) = G_TE(R) + 2*C0*real(QkE(I))*2*dk;
                G_TM(R) = G_TM(R) + 2*C0*real(QkM(I))*2*dk;
                R_sp(R) = R_sp(R) + 2*C0*real(Qksp(I))*2*dk;
                delta_n(R) = delta_n(R) - Consts.c*C0/omega(R)*imag(QkE(I))*2*dk;
            end
        end
    end
end

%%% -------------------------------------------------------------------------------------------------------------------------------%%%

function deltaE_sx = FindDeltaEsx(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, A_B, Ef, G, q_times_epsilon_q, dq, Epsr, T)

global Consts;

UT = Consts.k_B*T;
Co = (Consts.e_0/pi)^2/8/Consts.eps_0/Epsr;

theta = zeros(k_point, k_prime_point);
deltaE_sx = zeros(1,k_point);
Fermi = zeros(1,k_prime_point);

phi_min = 0;
phi_max = pi;
phi_point = 51;
dphi = (phi_max-phi_min)/(phi_point-1);

for (I=1:k_point)
    k = k_min+(I-1)*dk;
    deltaE_sx(I) = 0;
    for (J=1:k_prime_point)
        k_prime = k_prime_min+(J-1)*dk_prime;
        if (J<I)
            theta(I,J) = theta(J,I);
        else
            theta(I,J) = 0;
            if (I==1)
                Fermi(J) = 1/(1+exp((A_B(J)-Ef)/UT));
            end
            for (K=1:phi_point)
                phi = phi_min+(K-1)*dphi;
                q_prime = sqrt(k.^2+k_prime.^2-2.*k.*k_prime.*cos(phi));
                v1 = q_prime/dq;
                v2 = floor(v1);
                v3 = v1-v2;
                if (v2+2<k_point)
                    F = ((1-v3)*G(v2+1)+v3*G(v2+2))/((1-v3)*q_times_epsilon_q(v2+1)+v3*q_times_epsilon_q(v2+2));
                else
                    F = 0;
                end
                if ((K==1)||(K==phi_point))
                    theta(I,J) = theta(I,J)+Co*F*(2*dphi)/2;
                else
                    theta(I,J) = theta(I,J)+Co*F*(2*dphi);
                end
            end
        end
        if ((J==1)||(J==k_prime_point))
            deltaE_sx(I) = deltaE_sx(I) - k_prime*Fermi(J)*theta(I,J)*(dk_prime/2);
        else
            deltaE_sx(I) = deltaE_sx(I) - k_prime*Fermi(J)*theta(I,J)*dk_prime;
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

function theta = FindTheta(k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, G, q_times_epsilon_q, dq, Epsr_static)

global Consts;

phi_min = 0;
phi_max = pi;
phi_point = 51;
dphi=(phi_max-phi_min)/(phi_point-1);
Co = (Consts.e_0/pi)^2/8/Consts.eps_0/Epsr_static;

for (I=1:k_point)
    k = k_min+(I-1)*dk;
    for (J=1:k_prime_point)
        k_prime = k_prime_min+(J-1)*dk_prime;
        if (J<I)
            theta(I,J) = theta(J,I);
        else
            theta(I,J) = 0;
            for (K=1:phi_point)
                phi = phi_min + (K-1)*dphi;
                q_prime = sqrt(k.^2 + k_prime.^2 - 2*k*k_prime*cos(phi));
                v1 = q_prime/dq;
                v2 = floor(v1);
                v3 = v1-v2;
                if (v2+2 < k_point)
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

function [QkE, QkM, Qksp, Q] = Enhancement(theta, k_min, dk, k_point, k_prime_min, dk_prime, k_prime_point, ...
                                        E_c, Efc, E_v, Efv, T, omega, ukE, ukM, deltaE_ch, deltaE_sx_e, ...
                                        deltaE_sx_h, uk0_TE, uk0_TM, gamma, broad_type)

global Consts;

UT = Consts.k_B*T;
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
            k_prime(J) = k_prime_min+(J-1)*dk_prime;
            omega_k_prime(J) = (E_c(J)+E_v(J)+deltaE_ch+deltaE_sx_e(J)+deltaE_sx_h(J))/Consts.hbar;
            w_k_prime(J) = 1/(1+exp((E_c(J)-Efc)/UT))+1/(1+exp((E_v(J)-Efv)/UT))-1;
            if (strcmp(broad_type,'Lorentzian_enhanced'))
                if (omega<=omega_min)
                    RE(J) = (1-exp((Consts.hbar*omega-(Efc-E_c(1)+Efv+deltaE_ch+deltaE_sx_e(1)+deltaE_sx_h(1)))/UT))*(uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RM(J) = (1-exp((Consts.hbar*omega-(Efc-E_c(1)+Efv+deltaE_ch+deltaE_sx_e(1)+deltaE_sx_h(1)))/UT))*(uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                elseif ((omega>omega_min) && (omega<=omega_max))
                    RE1 = (1-exp((Consts.hbar*omega-(Efc-E_c(1)+Efv+deltaE_ch+deltaE_sx_e(1)+deltaE_sx_h(1)))/UT))*(uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RM1 = (1-exp((Consts.hbar*omega-(Efc-E_c(1)+Efv+deltaE_ch+deltaE_sx_e(1)+deltaE_sx_h(1)))/UT))*(uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
                    RE2 = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                    RM2 = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
                    RE(J) = RE1+(RE2-RE1)*(omega-omega_min)/(omega_max-omega_min);
                    RM(J) = RM1+(RM2-RM1)*(omega-omega_min)/(omega_max-omega_min);
                else
                    RE(J) = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                    RM(J) = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
                end
            else
                RE(J) = (uk0_TE*ukE(J)).^2*k_prime(J)*w_k_prime(J);
                RM(J) = (uk0_TM*ukM(J)).^2*k_prime(J)*w_k_prime(J);
            end
            Rsp(J) = 2*(uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT))/3 + ...
                       (uk0_TM*ukM(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT))/3;
            %Rsp(J) = (uk0_TE*ukE(J)).^2*k_prime(J)/(1+exp((E_c(J)-Efc)/UT))/(1+exp((E_v(J)-Efv)/UT));
        end
        
        %M = (theta(I,J)*k_prime(J)*w_k_prime(J)*dk_prime/Consts.hbar)/((omega_k_prime(J)-omega).^2+gamma^2);
        M = (theta(I,J)*w_k_prime(J)*dk_prime/Consts.hbar)/((omega_k_prime(J)-omega).^2+gamma^2);
        %M_full = (theta(I,J)*w_k_prime(J)*dk_prime/Consts.hbar)/((omega_k_prime(J)-omega).^2+gamma^2);
        Re_M_TE = Re_M_TE+M*ukE(J)/ukE(I)*(omega_k_prime(J)-omega);
        Im_M_TE = Im_M_TE+M*ukE(J)/ukE(I)*gamma;
        Re_M_TM = Re_M_TM+M*ukM(J)/ukM(I)*(omega_k_prime(J)-omega);
        Im_M_TM = Im_M_TM+M*ukM(J)/ukM(I)*gamma;
        %M_mat_TE(I,J) = (M_full).*(ukE(J)/ukE(I)*(omega_k_prime(J)-omega) + 1i*ukE(J)/ukE(I)*gamma);
        %M_mat_TM(I,J) = (M_full).*(ukM(J)/ukM(I)*(omega_k_prime(J)-omega) + 1i*ukM(J)/ukM(I)*gamma);
    end
    
    switch (broad_type)
        case 'Lorentzian',
            lineshape = gamma^2/((omega_k_prime(I)-omega)^2+gamma^2);
        case 'Lorentzian_enhanced',
            lineshape = gamma^2/((omega_k_prime(I)-omega)^2+gamma^2);
        case 'Sech',
            lineshape = 1/cosh((omega_k_prime(I)-omega)/gamma);
        case 'Cos_Hyp'
            lineshape = 1/cosh((omega_k_prime(I)-omega)/gamma);
    end
    Q = 1/(1+Re_M_TE+1i*Im_M_TE);
    QkE(I) = (RE(I)*lineshape/gamma^2)/(1+Re_M_TE+1i*Im_M_TE)*(gamma-1i*(omega_k_prime(I)-omega));
    QkM(I) = (RM(I)*lineshape/gamma^2)/(1+Re_M_TM+1i*Im_M_TM)*(gamma-1i*(omega_k_prime(I)-omega));
    Qksp(I) = (Rsp(I)*lineshape/gamma^2)/(1+2*Re_M_TE/3+Re_M_TM/3+1i*(2*Im_M_TE/3+Im_M_TM/3))*(gamma-1i*(omega_k_prime(I)-omega));
    %Qksp(I) = (Rsp(I)*lineshape/gamma^2)/(1+Re_M_TE+1i*(Im_M_TE))*(gamma-1i*(omega_k_prime(I)-omega));
end
