function [Q_k,M] = CalculateHFEnhancement(k_vec,omega, mu, w_tag, w, gamma, Theta, type)

global Consts;

k_len = length(k_vec);
Xi_0 = (-1i/Consts.hbar).*((omega.*mu)./(1i*(w_tag-w)+gamma));
if (strcmp(type,'Pade'))
	M = repmat((1./mu).', 1, k_len).*repmat(Xi_0.*k_vec, k_len, 1).*Theta;
    Q_k = 1./(1+trapz(k_vec, M, 2));
elseif (strcmp(type,'Exact'))
    d_k = k_vec(2)-k_vec(1);
    M = eye(k_len) - repmat((d_k./mu).', 1, k_len).*repmat(Xi_0, k_len, 1).*Theta;
    Q_k = M\ones(k_len,1);
elseif (strcmp(type,'Iter'))
    tol = 1e-3;
    Q_k = 0.01.*zeros(k_len,1);
    Q_temp = Q_k;
    while (1)
        Q_temp = Q_k;
        for (ii=1:k_len)
            Q_k(ii) = 1+(1/(mu(ii))).*trapz(k_vec, Theta(ii,:).*Xi_0.*Q_temp.');     
        end
        if (abs(sum(real(Q_k-Q_temp))) < 1e-4 && abs(sum(imag(Q_k-Q_temp))) < 1e-4)
           break; 
        end
    end
    M = 0;
end

%Q_k = abs(real(Q_k)) + 1i.*abs(imag(Q_k));
      