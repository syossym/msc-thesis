%
% This function calculates the screened dielectric function
%
function eps_q = Eps_q(k_vec,theta,E_fv,E_fc,T,f_e,f_h,E_e,E_h,eps_r)

global Consts;

%E_e = E_e - min(E_e);
%E_h = E_h - min(E_h);

%k_vec = k_vec(2:end);
k_len = length(k_vec);
% E_e = E_e(2:end);
% E_h = E_h(2:end);
% f_e = f_e(2:end);
% f_h = f_h(2:end);

for (kk = 1:k_len)
    int3 = zeros(size(k_vec)); int4 = int3;
    for (kk_tag = 1:k_len)
        k_tag = sqrt(k_vec(kk_tag)^2+k_vec(kk)^2-2*k_vec(kk_tag)*k_vec(kk)*cos(theta));
        
        k_tag(k_tag > max(k_vec)) = max(k_vec);
        k_tag(k_tag < min(k_vec)) = min(k_vec);
        
        E_e_k_minus_q(kk_tag, :) = interp1(k_vec, E_e, k_tag,'pchip');
        E_h_k_minus_q(kk_tag, :) = interp1(k_vec, E_h, k_tag,'pchip');
        E_h_k_minus_q(kk_tag, isnan(E_h_k_minus_q(kk_tag, :))) = E_h(kk);
        E_e_k_minus_q(kk_tag, isnan(E_e_k_minus_q(kk_tag, :))) = E_e(kk);
        f_h_k_minus_q(kk_tag, :) = FermiDirac(E_fv, E_h_k_minus_q(kk_tag, :), T);
        f_e_k_minus_q(kk_tag, :) = FermiDirac(E_fc, E_e_k_minus_q(kk_tag, :), T);
                
        nom_e(kk_tag, :) = f_e_k_minus_q(kk_tag, :) - f_e(kk_tag).*ones(size(theta));
        nom_h(kk_tag, :) = f_h_k_minus_q(kk_tag, :) - f_h(kk_tag).*ones(size(theta));
        den_e(kk_tag, :) = E_e_k_minus_q(kk_tag, :) - E_e(kk_tag).*ones(size(theta));
        den_h(kk_tag, :) = E_h_k_minus_q(kk_tag, :) - E_h(kk_tag).*ones(size(theta));
               
        int1(kk_tag, :) = nom_e(kk_tag, :)./den_e(kk_tag, :);
        int2(kk_tag, :) = nom_h(kk_tag, :)./den_h(kk_tag, :);     
        %int1(kk_tag, isnan(int1(kk_tag,:))) = min(int1(kk_tag,:)); 
        %int2(kk_tag, isnan(int2(kk_tag,:))) = min(int2(kk_tag,:));

%         figure(1);
%         subplot(3,2,[1 2]);  plot(k_tag, E_e_k_minus_q(kk_tag, :), 'r.', k_vec, E_e, 'b-',...
%                             k_tag, E_h_k_minus_q(kk_tag, :), 'r.', k_vec, E_h, 'g-');
%         subplot(3,2,[3,4]);  plot(k_tag, f_e_k_minus_q(kk_tag, :), 'r.', k_vec, f_e, 'b-',...
%                             k_tag, f_h_k_minus_q(kk_tag, :), 'r.', k_vec, f_h, 'g-');
%         subplot(325); hold on;  plot(theta, int1(kk_tag, :), 'b');
%         subplot(326); hold on;  plot(theta, int2(kk_tag,:), 'r');
%         drawnow;  

        %p1 = polyfit(theta, int1(kk_tag,:), 25);
        %p2 = polyfit(theta, int2(kk_tag,:), 25);
        int3(kk_tag) = trapz(theta, (int1(kk_tag, :)));  %polyval(p1,theta));
        int4(kk_tag) = trapz(theta, (int2(kk_tag, :)));  %polyval(p2,theta));
    end
    
%     int1(1,:) = int1(2,:);
%     int2(1,:) = int2(2,:);
%     int1(k_len,:) = int1(k_len-1,:);
%     int2(k_len,:) = int2(k_len-1,:);
%     
%     for (kk_tag = 1:k_len)
%         int3(kk_tag) = trapz(theta, (int1(kk_tag, :)));  %polyval(p1,theta));
%         int4(kk_tag) = trapz(theta, (int2(kk_tag, :)));  %polyval(p2,theta));
%     end
    
    int3(isnan(int3))=0; int4(isnan(int4))=0;
    
    %p3 = polyfit(k_vec, int3, 20);
    %p4 = polyfit(k_vec, int4, 20);
    
    int3(end) = int3(end-1); int4(end) = int4(end-1);
    
%     figure(5);
%     subplot(211); plot(k_vec, int3, 'r.', k_vec, smooth(int3, 'rloess')); %polyval(p3,k_vec), 'b');
%     title([num2str(kk)]);
%     subplot(212); plot(k_vec, int4, 'r.', k_vec, smooth(int4, 'rloess')); %polyval(p4,k_vec), 'b');
%     drawnow;
    
    %int3 = smooth(int3, 'rloess'); int4 = smooth(int4, 'rloess');  %polyval(p3,k_vec);  int4 =  polyval(p4,k_vec);
    alpha_q(kk) = (trapz(k_vec, k_vec.*int3) + trapz(k_vec,  k_vec.*int4));
    
end
% for (ii=1:k_len)
%     if (isnan(alpha_q(ii)))
%         if (ii-1<1)
%             alpha_q(ii) = alpha_q(ii+1)/2;
%         elseif (ii+1>k_len)
%             alpha_q(ii) = 1.5*alpha_q(ii-1);
%         else
%             alpha_q(ii) = alpha_q(ii-1) + (alpha_q(ii+1)-alpha_q(ii-1))/2;
%         end
%     end
% end
% alpha_q = [0, alpha_q];
%p_alpha = smooth(p_alpha); %polyfit(k_vec, alpha_q, 20);
%eps_q = k_vec - (Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)).*alpha_q; %polyval(p_alpha,k_vec);
eps_q = (Consts.e_0^2./(4*pi^2*eps_r*Consts.eps_0)).*alpha_q;