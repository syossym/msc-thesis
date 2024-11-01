function E_f = QuasiFermiLevels(V_max, subband_min, DOS, N, T)

%
% Solutions are sought, using a stepping algorithm, a midpoint rule and 
% finally a Newton-Raphson method, to the equation f(E_F)=0.  This 
% function has been derived by integrating the total density of states
% multiplied by the Fermi-Dirac distribution function across the in-plane
% energy thus giving the total charge carrier density, i.e.
% 
%      oo
% Ne=  I  f(E)N(E) dE
%      Eo
% 
% where f(E) is the normal Fermi-Dirac distribution function and
% N(E)=m/(pi*sqr(hbar)) is the areal density of states in a QW.
%
%   Input: 'V_max' - the potential top of QW
%          'subband_min' (J) - the reference energy for the calculation
%          'DOS' - the density of states for the calculation
%          'N' - the total charge carrier density
%          'T' - the simualation temperature
%           
%   Output: 'E_f' (J) - the calculated quasi Fermi level energy
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, October 2009
% Edited by: -
%

global Consts;

dE = 0.01*1e-3*Consts.e_0;

E_min = subband_min;
%V_max = max(V);

x = E_min - 20*Consts.k_B*T;    % first value of x

% In this implementation, the upper limit of integration is set at the 
% Fermi level+10kT, limited at potential maximum.

E_max = x + 10*Consts.k_B*T;      % subband maximum (top of QW)
if (E_max > V_max)
    E_max = V_max;
end

y2 = f(x,E_max,E_min,DOS,N,T);

while(1)
    y1 = y2;
    x = x + dE;
    E_max = x + 10*Consts.k_B*T;     
    if (E_max > V_max)
        E_max = V_max;
    end
    y2 = f(x,E_max,E_min,DOS,N,T);

    if (y1*y2<=0)
        break;
    end
end

% Improve estimate using midpoint rule
x = x - abs(y2)/(abs(y1)+abs(y2))*dE;

E_f = x;

% -------- Inner Functions --------- %

function y = f(E_f,E_max,E_min,DOS,N,T)

global Consts; 

m = (DOS*pi*Consts.hbar^2);
y = ((m*Consts.k_B*T)/(pi*Consts.hbar^2))*...
		(...
		((E_max-E_f)/(Consts.k_B*T)-log(1+exp((E_max-E_f)/(Consts.k_B*T))))-...
		((E_min-E_f)/(Consts.k_B*T)-log(1+exp((E_min-E_f)/(Consts.k_B*T))))...
		)...
	 -N;
    