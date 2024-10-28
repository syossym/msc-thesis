function E_f = QuasiFermiLevels(sol_method,V,k_t,E_k,N,T,L_z,DOS)

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

switch (sol_method)
    
    case 'Iterative',
        
        dE = 0.01*1e-3*Consts.e_0;
        
        E_min = min(E_k);
        V_max = max(V);
        
        x = E_min - 20*Consts.k_B*T;    % first value of x
        
        % In this implementation, the upper limit of integration is set at the
        % Fermi level+10kT, limited at potential maximum.
        
        E_max = x + 10*Consts.k_B*T;      % subband maximum (top of QW)
        if (E_max > V_max)
            E_max = V_max;
        end
        
        y2 = target_f_1(x,E_max,E_min,DOS,N,T);
        
        while(1)
            y1 = y2;
            x = x + dE;
            E_max = x + 10*Consts.k_B*T;
            if (E_max > V_max)
                E_max = V_max;
            end
            y2 = target_f_1(x,E_max,E_min,DOS,N,T)
            
            if (y1*y2<=0)
                break;
            end
        end
        
        % Improve estimate using midpoint rule        x = x - abs(y2)/(abs(y1)+abs(y2))*dE;
        
        E_f = x;
        
    case 'fsolve',
        
        options=optimset('TolFun',1e-100,'TolX',1e-100,'MaxIter',5000,'MaxFunEvals',5000,'Display','iter');

         min(E_k{6})
        
        %x = fsolve(@(E_f) target_f_2(k_t, E_k, E_f, N, T, L_z), min(E_k{1}), options);
        
        %res = target_f_2(k_t, E_k, x, N, T, L_z)
        E_f = 0; r_vec = [];
        while(E_f < max(V))
           E_f;
           res = target_f_2(k_t, E_k, E_f, N, T, L_z);
           r_vec = [r_vec,res];
           
           if (abs(res) < 1e-1)
              break; 
           end
           
            E_f = E_f + 0.00001*Consts.e_0;
        end
                
%         if (res ~= 0)
%             error('QuasiFermiLevels:convergence_error', 'Error');
%         end
        E_f
        %E_f = x;
        
    otherwise
       
        error('QuasiFermiLevels:input_argument_error',...
            'The solution method not defined correctly.');
        
end

% -------- Inner Functions --------- %

function y = target_f_1(E_f,E_max,E_min,DOS,N,T)

global Consts;

m = (DOS*pi*Consts.hbar^2);
y = ((m*Consts.k_B*T)/(pi*Consts.hbar^2))*...
    (...
    ((E_max-E_f)/(Consts.k_B*T)-log(1+exp((E_max-E_f)/(Consts.k_B*T))))-...
    ((E_min-E_f)/(Consts.k_B*T)-log(1+exp((E_min-E_f)/(Consts.k_B*T))))...
    )...
    -N;

function y = target_f_2(k_t, E_k, E_f, N, T, L_z)

global Consts;

sum = 0;
for (ii=1:length(E_k))
    sum = sum + trapz(k_t, k_t./(exp((E_k{ii}-E_f)/(Consts.k_B*T)) + 1));
end
(1/(pi*L_z)).*sum;
y = (1/(pi*L_z)).*sum - N;


