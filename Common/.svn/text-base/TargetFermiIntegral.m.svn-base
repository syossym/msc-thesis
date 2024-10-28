
function y = TargetFermiIntegral(k_t, E_k, E_f, N, T, L_z, type)

global Consts;

sum = 0;
for (ii=1:length(E_k))
    if (strcmp(type,'c'))
        I = trapz(k_t, k_t./(exp((E_k{ii}-E_f)/(Consts.k_B*T)) + 1));
    elseif (strcmp(type,'v'))
        if (E_k{1}(1)>0)
            for (ii=1:length(E_k))
                E_k{ii} = -E_k{ii};
            end
        end
        I = 2.*trapz(k_t, k_t.*(1-(exp((E_k{ii}-E_f)/(Consts.k_B*T)) + 1).^(-1)));
    elseif (strcmp(type,'h'))
        if (E_k{1}(1)<0)
            for (ii=1:length(E_k))
                E_k{ii} = -E_k{ii};
            end
        end
        I = 2*trapz(k_t, k_t./(exp((E_k{ii}-E_f)/(Consts.k_B*T)) + 1));
    end
    sum = sum + I;
end

y = N/L_z - (1/(pi*L_z)).*sum;
