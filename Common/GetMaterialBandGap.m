function E_g = GetMaterialBandGap(x,T)

switch (x)
    case 0,  % GaAs
        E_g = 1.521 - 5.405e-4*T^2/(T+204);    % Gamma [eV]
    case 1,  % AlAs
        E_g = 2.239 - 6e-4*T^2/(T+408);        % X [eV]
    otherwise,  % GaAlAs
        E_g_GaAs_X = 1.981-(4.6e-4*T^2)/(204+T);
        E_g_AlAs_X = 2.239-(6e-4*T^2)/(408+T);
        E_g_GaAs_G = 1.521-(5.58e-4*T^2)/(220+T);
        E_g_AlAs_G = 2.891-(8.78e-4*T^2)/(332+T);
        
        E_g_X = (1-x)*E_g_GaAs_X + x*E_g_AlAs_X - 0.143*(1-x)*x;
        E_g_G = (1-x)*E_g_GaAs_G + x*E_g_AlAs_G;
        
        E_g = min([E_g_X, E_g_G]);    % [eV]
end