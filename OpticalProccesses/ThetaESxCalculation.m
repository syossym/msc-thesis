function [Theta, Theta_1, E_SX] = ThetaESxCalculation(Gq_ee,Gq_hh,Gq_eh,theta,eps_q,k,f_e,f_h,eps_r)

global Consts;

C = Consts.e_0^2./(8*pi^2*eps_r*Consts.eps_0);
E_SX = zeros(1,length(k)); int_E_SX = zeros(1,length(k));
Theta = zeros(length(k));
Theta_1 = zeros(length(k));
for (kk = 1:length(k))
    for (kk_tag = 1:length(k))
        q_tag = sqrt(k(kk).^2 + k(kk_tag).^2 - 2.* k(kk).*k(kk_tag).*cos(theta));
        
        q_tag(q_tag > max(k)) = max(k);
        q_tag(q_tag < min(k)) = min(k);
        
        % Theta
        form_eh = interp1(k,Gq_eh,q_tag,'pchip');
        if (sum(isnan(form_eh)) ~= 0)
            form_eh(isnan(form_eh)) = Gq_eh(kk);
        end
        eps = interp1(k,eps_q,q_tag,'pchip');
        if (sum(isnan(eps)) ~= 0)
            eps(isnan(eps)) = eps_q(kk);
        end
        Theta(kk,kk_tag) = C.*trapz(theta, form_eh./eps);
        Theta_1(kk,kk_tag) = C.*trapz(theta, form_eh./q_tag);
        
        % E_SX
        form_ee = interp1(k,Gq_ee,q_tag,'pchip');
        if (sum(isnan(form_ee)) ~= 0)
            form_ee(isnan(form_ee)) = Gq_ee(kk);
        end
        %hold on; plot(theta, form_ee, 'r'); hold off;
        form_hh = interp1(k,Gq_hh,q_tag,'pchip');
        if (sum(isnan(form_hh)) ~= 0)
            form_hh(isnan(form_hh)) = Gq_hh(kk);
        end
        %hold on; plot(theta, form_hh, 'r'); hold off;
        integ = (form_ee.*f_e(kk_tag)+form_hh.*f_h(kk_tag))./eps;
        %p1 = polyfit(theta, integ, 5);
        
%         if (kk>60)
%             figure(3);
%             subplot(221); plot(k, Gq_ee, 'b', q_tag, form_ee, '.r'); title('G_q^{ee}');
%             subplot(222); plot(k, Gq_hh, 'b', q_tag, form_hh, '.r'); title('G_q^{hh}');
%             subplot(223); plot(k, eps_q, 'b', q_tag, eps, '.r'); title('\epsilon_q');
%             subplot(224); plot(q_tag, integ);
%         end
        
        int_E_SX(kk_tag) = trapz(theta, integ); %polyval(p1,theta));
    end
    %p2 = polyfit(k, int_E_SX, 15);
    %int_E_SX = RepairNanValues(int_E_SX); %polyval(p2,k));
    %int_E_SX = smooth(k,int_E_SX,0.8,'rloess').';
    %Theta(kk,:) = smooth(k,Theta(kk,:),0.8,'rloess').';
   
%     figure(4);
%     subplot(211);
%     plot(k, int_E_SX, 'r'); title([num2str(kk)]); 
%     subplot(212);
%     plot(k, Theta(kk,:), 'r');
%     drawnow;
    
    E_SX(kk) = -C.*trapz(k, k.*int_E_SX);
end
% p3 = polyfit(k,E_SX, 15);
% E_SX = polyval(p3,k);