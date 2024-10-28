function n = ModifiedSellmeierEquation(x, T, E)
%
% This function calculates the refractive index of a Al{x}Ga{1-x}As 
% alloy according to the empirical Sellmeier equation, as presented in
% OPTICS LETTERS, Vol. 32, No. 5 (March 1, 2007) and the website
% http://www.batop.com/information/n_GaAs.html.
%
%   Input: 'x' - Al molar fraction in the alloy,
%          'T' - the ambience temperature (K),
%          'E' - the energy vector for the calculation (J).
%
%   Output: 'n' - the result refractive index.
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, May 2010
% Edited by: -
%

global Consts;

if (x <= 0.36)
   C = (0.52886-0.735*x)^2; 
else
   C = (0.30386-0.105*x)^2;
end

lambda = 1e6*(2*pi./E)*Consts.c*Consts.hbar;    % [m]
T_C = T - 273.15;                           % [oC]

n = sqrt(10.906-2.92*x+(0.97501./(lambda.^2+C))-0.002467*(1.41*x+1).*lambda.^2)+ ...
    ((T-26)*(2.04-0.3*x)*1e-4);

% lambda = (2*pi./E)*Consts.c*Consts.hbar;    % [m]
% 
% if (x == 0)              % GaAs 
%     A = 8.950;
%     B = 2.054;
%     C_2 = 0.390;
%     
%     n = sqrt(A+B./(1-C_2./(lambda*1e6).^2));   
% elseif (x == 1)          % AlAs
%     A0 = 25.3;
%     B0 = -0.8;
%     E0 = 2.95*Consts.e_0;             % [J]
%     E0_p_Delta0 = 3.25*Consts.e_0;    % [J]
%     
%     chi = 2*pi*Consts.hbar*Consts.c./(lambda*E0);
%     chi_SO = 2*pi*Consts.hbar*Consts.c./(lambda*E0_p_Delta0);
%     f = (2-sqrt(1+chi)-sqrt(1-chi))./chi.^2;
%     f_SO = (2-sqrt(1+chi_SO)-sqrt(1-chi_SO))./chi_SO.^2;
%     n = sqrt(A0*(f+(f_SO./2).*(E0./E0_p_Delta0)^(3/2))+B0);
% else                     % AlGaAs
%     A0 = 6.3 + 19.0*x;
%     B0 = 9.4 - 10.2*x;
%     E0 = (1.425 + 1.155*x + 0.37*x^2)*Consts.e_0;             % [J]
%     E0_p_Delta0 = (1.765 + 1.115*x + 0.37*x^2)*Consts.e_0;    % [J]
%     
%     chi = E./E0;
%     chi_SO = E./E0_p_Delta0;
%     f = (2-sqrt(1+chi)-sqrt(1-chi))./chi.^2;
%     f_SO = (2-sqrt(1+chi_SO)-sqrt(1-chi_SO))./chi_SO.^2;
%     n = sqrt(A0*(f+(f_SO./2).*(E0./E0_p_Delta0)^(3/2))+B0);
% end
