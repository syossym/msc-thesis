function n = CalculateRefractiveIndex(method, x, T, E)

global Consts aquila_control

if (x==1)
    material = 'AlAs';
elseif (x==0)
    material = 'GaAs';
else
    material = 'GaAlAs';
end

switch (method)
    case 'exp_1',
        n = GetExpRefractiveIndex(x, T, E);
    case 'exp_2',
        n = GetRefractiveIndex(material, E, x);
    case 'emp_1',
        n = EmpiricalRefraciveIndex('Sellemeier',x,T,E);
    case 'emp_2',
        n = EmpiricalRefraciveIndex('Zhang',x,T,E);
end

function n = GetExpRefractiveIndex(x, T, E_interp)

global Consts;

% Load the experimental data
n_GaAs = load('Materials\n_GaAs_300K.mat');
n_AlAs = load('Materials\n_AlAs_300K.mat');
n_GaAlAs = load('Materials\n_GaAlAs_300K.mat');

% Calculate the refractive index
if (x==1)
    d_Eg = GetMaterialBandGap(x,T)-GetMaterialBandGap(x,300);
    n_real = interp1(n_AlAs.E+d_Eg, n_AlAs.AlAs(:,1).', E_interp, 'pchip');
    n_imag = interp1(n_AlAs.E+d_Eg, n_AlAs.AlAs(:,2).', E_interp, 'pchip');
    %n_imag(n_AlAs.E+d_Eg<GetMaterialBandGap(x,T)) = 0;
    n = n_real + 1i*n_imag;
elseif (x==0)
    d_Eg = GetMaterialBandGap(x,T)-GetMaterialBandGap(x,300);;
    n_real = interp1(n_GaAs.E+d_Eg, n_GaAs.GaAs(:,1).', E_interp, 'pchip');
    n_imag = interp1(n_GaAs.E+d_Eg, n_GaAs.GaAs(:,2).', E_interp, 'pchip');
    
    %n_imag(n_imag<0.05) = 0;
    n = n_real + 1i*n_imag;
else
    %     if (x < 0.69)
    %         aquila_control.T = T;
    %         E_g_T = gaasmaterial(x,'E_G6G8');
    %         aquila_control.T = 300;
    %         E_g_300 = gaasmaterial(x,'E_G6G8');
    %         d_Eg = E_g_T-E_g_300;
    %     else
    
    %d_Eg = GetMaterialBandGap(x,T)-GetMaterialBandGap(x,300);
    %end
    
    %     E = n_AlAs.E;
    %     n_real = (1-x).*(n_GaAs.GaAs(1:length(E),1)) + x.*(n_AlAs.AlAs(:,1));
    %     n_imag = (1-x).*(n_GaAs.GaAs(1:length(E),2)) + x.*(n_AlAs.AlAs(:,2));
    %
    %     %n_imag(n_imag<0.05) = 0;
    %
    %     n_real = interp1(E+d_Eg, n_real, E_interp, 'pchip');
    %     n_imag = interp1(E+d_Eg, n_imag, E_interp, 'pchip');
    %     n = n_real + 1i*n_imag;
    
    d_Eg = GetMaterialBandGap(x,T)-GetMaterialBandGap(x,300);
    n_300K = n_GaAlAs.AlGaAs{floor(x*10)};
    E = 1:n_300K.E(2)-n_300K.E(1):1.5-(n_300K.E(2)-n_300K.E(1));
    n_imag = [zeros(1,length(E)), n_300K.n_imag.'];
    E = [E, n_300K.E];
    n_real = interp1(n_300K.E+d_Eg, n_300K.n_real.', E_interp, 'cubic');
    n_imag = interp1(E+d_Eg, n_imag, E_interp, 'pchip');
    n_imag(E_interp<1.5) = 0;
    n = n_real + 1i*n_imag;
    
end

function n = EmpiricalRefraciveIndex(method,x,T,E)

global Consts;

switch (method)
    case 'Sellemeier',
        n = ModifiedSellmeierEquation(x, T, E);
    case 'Zhang',
        switch (x)
            case 0,    % GaAs
                E_0 = 1.5192 + 2.862e-2*(1-coth(1.59e-2/(2*Consts.k_B_eV*T))) + ...
                    3.696e-2*(1-coth(3.36e-2/(2*Consts.k_B_eV*T)));
                A = 7.3377 + 5.534e-4*T - 3.56e-7*T^2;
                E_1_sq = 3.791 - 3.779e-4*T - 1.121e-6*T^2;
                
                n = sqrt(A+0.001680081./(E_0^2-E.^2)+13.603615./(E_1_sq-E.^(1.22)));
            case 1,    % AlAs
                E_0 = 3.099 + 2.862e-2*(1-coth(1.59e-2/(2*Consts.k_B_eV*T))) + ...
                    3.696e-2*(1-coth(3.36e-2/(2*Consts.k_B_eV*T)));
                A = 2.857 + 4.574e-4*T - 2.942e-7*T^2;
                E_1_sq = 11.717 - 3.779e-4*T - 1.121e-6*T^2;
                R = 2.157e-3./(1.331e-3-E.^2);
                
                n = sqrt(A + 0.060876./(E_0^2-E.^2) + 61.064215./(E_1_sq-E.^2)+R);
            otherwise, % GaAlAs
                n_GaAs = EmpiricalRefraciveIndex('Zhang',0,T,E);
                n_AlAs = EmpiricalRefraciveIndex('Zhang',1,T,E);
                n = (1-x).*n_GaAs + x.*n_AlAs;
        end
end