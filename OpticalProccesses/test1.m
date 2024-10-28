L_z = 200e-10;
T = 300;
gamma = 1e13;
init_E = 0.01*Consts.e_0;
dk = R.k_t_mat(2) - R.k_t_mat(1);
E_exc = (R.Materials{2}.E_g-0.1)*Consts.e_0 : 1e-3*Consts.e_0 : (R.Materials{2}.E_g+0.5)*Consts.e_0;
Mb_s = Consts.e_0*R.Materials{2}.E_p*Consts.m_0/6;
C_g = 1./(Consts.hbar*Consts.c*Consts.eps_0*L_z*pi*R.Materials{2}.n);
k_vec = R.k_t_mat(2:end);
k_len = length(k_vec);
k_mat = repmat(k_vec.^2, k_len, 1);
theta = [0:0.1:2*pi];
z_exp = -abs(repmat(R.z_grid,1,length(R.z_grid)) -...
    repmat(R.z_grid',length(R.z_grid),1));

N = 2e11;
P = N
cb_num = 1; vb_num = 1;

% Quasi Fermi levels
options = optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');
E_fv = fzero(@(E_fv) TargetFermiIntegral(R.k_t_mat, R.E_k_v, E_fv, P*1e4, T, L_z), init_E, options);
E_fc = fzero(@(E_fc) TargetFermiIntegral(R.k_t_mat, R.E_k_c, E_fc, N*1e4, T, L_z), init_E, options);

E_e = R.E_k_c{cb_num}(1:k_len);
E_h = R.E_k_v{vb_num}(1:k_len);
f_h = FermiDirac(E_fv, E_h, T);
f_e = FermiDirac(E_fc, E_e, T);
mu_TE = sqrt(squeeze(R.rtsTE(vb_num,cb_num,(1:k_len)))).'.*Mb_s;
mu_TM = sqrt(squeeze(R.rtsTM(vb_num,cb_num,(1:k_len)))).'.*Mb_s;
Psi_e = abs(R.wf_e_mat{cb_num}(:,2));
Psi_hh = abs(R.wf_h_mat{vb_num}(:,2));
Psi_lh = abs(R.wf_h_mat{vb_num}(:,3));
Xi_e = repmat(Psi_e.^2, 1, length(R.z_grid));
Xi_h = repmat(Psi_hh.'.^2 + Psi_lh.'.^2, length(R.z_grid), 1);
Xi_ee = Xi_e.*Xi_e.';
Xi_hh = Xi_h.'.*Xi_h;
Xi_eh = Xi_e.*Xi_h;

k_vec_interp = 0:1*dk:k_vec(end);
f_e_interp = interp1(k_vec, f_e, k_vec_interp);
f_h_interp = interp1(k_vec, f_h, k_vec_interp);
E_e_interp = interp1(k_vec, E_e, k_vec_interp);
E_h_interp = interp1(k_vec, E_h, k_vec_interp);

for (kk=1:k_len)
    kk
    % Gq
    z_exp_q = exp(z_exp.*k_vec(kk));
    G_ee(kk) = trapz(R.z_grid, trapz(R.z_grid,Xi_ee.*z_exp_q,1), 2);
    G_hh(kk) = trapz(R.z_grid, trapz(R.z_grid,Xi_hh.*z_exp_q,1), 2);
    G_eh(kk) = trapz(R.z_grid, trapz(R.z_grid,Xi_eh.*z_exp_q,1), 2);
    
    % q*eps_q
    int3 = zeros(size(k_vec)); int4 = int3;
    for (kk_tag = 1:k_len)
        kk_tag;
        k_tag(kk_tag,:) = sqrt(k_vec(kk_tag)^2+k_vec(kk)^2-2*k_vec(kk_tag)*k_vec(kk)*cos(theta));
        
        E_e_k_minus_q(kk_tag, :) = interp1(k_vec_interp, E_e_interp, k_tag(kk_tag, :));
        
%         E_e_k_minus_q(kk_tag, :) = R.Materials{2}.E_g*Consts.e_0 + ...
%                                    R.E_0_c(cb_num).*1e-3*Consts.e_0 + ...
%                                    (Consts.hbar^2*k_tag(kk_tag,:).^2)./(2*R.Materials{2}.m_e);
        E_h_k_minus_q(kk_tag, :) = interp1(k_vec_interp, E_h_interp, k_tag(kk_tag, :));
        
        E_h_k_minus_q(kk_tag, isnan(E_h_k_minus_q(kk_tag, :))) = E_h(kk);
        E_e_k_minus_q(kk_tag, isnan(E_e_k_minus_q(kk_tag, :))) = E_e(kk);
        f_h_k_minus_q(kk_tag, :) = FermiDirac(E_fv, E_h_k_minus_q(kk_tag, :), T);
        f_e_k_minus_q(kk_tag, :) = FermiDirac(E_fc, E_e_k_minus_q(kk_tag, :), T);
%         f_e_k_minus_q(kk_tag, :) = interp1(k_vec_interp, f_e_interp, k_tag(kk_tag, :));
%         f_h_k_minus_q(kk_tag, :) = interp1(k_vec_interp, f_h_interp, k_tag(kk_tag, :));
        
   
%         f_e_k_minus_q(kk_tag, isnan(f_e_k_minus_q(kk_tag, :))) = f_e(kk);
%         f_h_k_minus_q(kk_tag, isnan(f_h_k_minus_q(kk_tag, :))) = f_h(kk);
%         E_e_k_minus_q(kk_tag, isnan(E_e_k_minus_q(kk_tag, :))) = E_e(kk);
%         E_h_k_minus_q(kk_tag, isnan(E_h_k_minus_q(kk_tag, :))) = E_h(kk);

        nom_e(kk_tag, :) = f_e_k_minus_q(kk_tag, :) - f_e(kk_tag).*ones(size(theta));
        nom_h(kk_tag, :) = f_h_k_minus_q(kk_tag, :) - f_h(kk_tag).*ones(size(theta));
        den_e(kk_tag, :) = E_e_k_minus_q(kk_tag, :) - E_e(kk_tag).*ones(size(theta));
        den_h(kk_tag, :) = E_h_k_minus_q(kk_tag, :) - E_h(kk_tag).*ones(size(theta));
        
%         figure(1); 
%         subplot(221); plot(k_tag(kk_tag,:), f_e(kk_tag).*ones(size(theta)), ':b', k_tag(kk_tag,:), f_e_k_minus_q(kk_tag, :), 'r');
%         subplot(222); plot(k_tag(kk_tag,:), f_h(kk_tag).*ones(size(theta)), ':b', k_tag(kk_tag,:), f_h_k_minus_q(kk_tag, :), 'r');        
%         subplot(223); plot(k_tag(kk_tag,:), E_e(kk_tag).*ones(size(theta)), ':b', k_tag(kk_tag,:), E_e_k_minus_q(kk_tag, :), 'r');        
%         subplot(224); plot(k_tag(kk_tag,:), E_h(kk_tag).*ones(size(theta)), ':b', k_tag(kk_tag,:), E_h_k_minus_q(kk_tag, :), 'r');
%         drawnow;
        
        int1(kk_tag, :) = nom_e(kk_tag, :)./den_e(kk_tag, :);
        int2(kk_tag, :) = nom_h(kk_tag, :)./den_h(kk_tag, :);
        
        int1(kk_tag, isnan(int1(kk_tag,:))) = 0; int2(kk_tag, isnan(int2(kk_tag,:))) = 0;
        
        p1 = polyfit(theta, int1(kk_tag,:), 25);
        p2 = polyfit(theta, int2(kk_tag,:), 25);
        
        int3(kk_tag) = k_vec(kk_tag).*trapz(theta, polyval(p1,theta));
        int4(kk_tag) = k_vec(kk_tag).*trapz(theta, polyval(p2,theta));
        
%         figure(1);
%         plot(k_tag(kk_tag,:), f_e_k_minus_q(kk_tag, :));
        
%         figure(2); title(['k=' num2str(k_vec(kk)) ', k''=' num2str(k_vec(kk_tag))]);
%         subplot(221); plot(k_vec, f_e, 'b', k_tag(kk_tag,:), f_e_k_minus_q(kk_tag,:), 'r.');
%         subplot(222); plot(k_vec, f_h, 'b', k_tag(kk_tag,:), f_h_k_minus_q(kk_tag,:), 'r.');
%         subplot(223); plot(k_vec, E_e, 'b', k_tag(kk_tag,:), E_e_k_minus_q(kk_tag,:), 'r.');
%         subplot(224); plot(k_vec, E_h, 'b', k_tag(kk_tag,:), E_h_k_minus_q(kk_tag,:), 'r.');            
%         drawnow;
%         figure(3);
%         subplot(221); plot(k_tag(kk_tag, :), nom_e(kk_tag, :), 'r.');
%         subplot(222); plot(k_tag(kk_tag, :), nom_h(kk_tag, :), 'r.');
%         subplot(223); plot(k_tag(kk_tag, :), den_e(kk_tag, :), 'r.');
%         subplot(224); plot(k_tag(kk_tag, :), den_h(kk_tag, :), 'r.');
%         drawnow;
%         figure(4);
%         subplot(211); plot(k_tag(kk_tag, :), int1(kk_tag, :), 'r.');
%         subplot(212); plot(k_tag(kk_tag, :), int2(kk_tag, :), 'r.');
%         drawnow;

    end
    
    int3(isnan(int3))=0; int4(isnan(int4))=0;
    
    p3 = polyfit(k_vec, int3, 20);
    p4 = polyfit(k_vec, int4, 20);
    
    figure(5);
    subplot(211); plot(k_vec, int3, 'r.', k_vec, polyval(p3,k_vec), 'b');
    subplot(212); plot(k_vec, int4, 'r.', k_vec, polyval(p4,k_vec), 'b');
    drawnow;
    
    int3 =  polyval(p3,k_vec);  int4 =  polyval(p4,k_vec);
    
    %         figure(1);
    %         
    %         drawnow;
    
%     figure(1);
%     subplot(221); plot(k_vec, f_e, 'bo', k_tag, f_e_k_minus_q, 'r.');
%     subplot(222); plot(k_vec, f_h, 'bo', k_tag, f_h_k_minus_q, 'r.');
%     subplot(223); plot(k_vec, E_e, 'bo', k_tag, E_e_k_minus_q, 'r.');
%     subplot(224); plot(k_vec, E_h, 'bo', k_tag, E_h_k_minus_q, 'r.');
%     drawnow;
%     
%     nom_e = f_e_k_minus_q - repmat(f_e.',1,length(theta));
%     den_e = E_e_k_minus_q - repmat(E_e.',1,length(theta));
%     nom_h = f_h_k_minus_q - repmat(f_h.',1,length(theta));
%     den_h = E_h_k_minus_q - repmat(E_h.',1,length(theta));
%     
%     nom_e(isnan(nom_e)) = 0; nom_h(isnan(nom_e)) = 0; 
%     den_e(isnan(nom_e)) = 0; den_h(isnan(nom_e)) = 0;
%     
%     figure(2);
%     subplot(221); plot(k_tag, nom_e, 'r');
%     subplot(222); plot(k_tag, nom_h, 'r');
%     subplot(223); plot(k_tag, den_e, 'r');
%     subplot(224); plot(k_tag, den_h, 'r');
%     drawnow;

%     int1 = (f_e_k_minus_q - repmat(f_e.',1,length(theta)))./(E_e_k_minus_q - repmat(E_e.',1,length(theta)));
%     int2 = (f_h_k_minus_q - repmat(f_h.',1,length(theta)))./(E_h_k_minus_q - repmat(E_h.',1,length(theta)));
%     int1(isnan(int1)) = 0; int2(isnan(int2)) = 0;
%     int3 = k_vec.'.*trapz(theta, int1, 2);
%     int4 = k_vec.'.*trapz(theta, int2, 2);
    
    %(Consts.e_0^2./(4*pi^2*R.Materials{2}.eps_r)).*(trapz(k_vec, int3.') + trapz(k_vec, int4.'))
    
    q_alpha_q(kk) = (Consts.e_0^2./(4*pi^2*R.Materials{2}.eps_r*Consts.eps_0)).*(trapz(k_vec, int3) + trapz(k_vec, int4));
end

p_alpha = polyfit(k_vec, q_alpha_q, 20);
q_eps_q = k_vec - polyval(p_alpha,k_vec);
