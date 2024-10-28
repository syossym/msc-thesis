function [a, E, error] = WavefunctionMinimum(E_min, E_max, k_t, E_acc, E_step, min_num)
% This function looks for a minimum in the wavefunction amplitude between
% E_min and E_max, at intervals of E_step, and for wavevector of k_t. If
% a minimum is found, its determined with an accuracy E_acc. 
% Optionally, the function can look for the min_num-th minimum,
% disregarding the first mun_num-1 minima in the interval. 
% The return values are a, which is the free parameter in the initial
% condition, the candidate eigenenergy E, and a variable error, which is
% set to 1 if no minimum is found and 0 otherwise.

% Global variables
global z  hbar  meV  V  dz  N  gamma1  gamma2  gamma3;

% Checking input
dE = E_max - E_min;

if ~exist('E_step')
    E_step = dE / 10;
end

if ~exist('min_num')
    min_num = 1;
end

if (E_step < E_acc)
    E_step = E_acc;
    fl = 1;
end

%% Finding the minimum

E           = E_min - 3*E_step;
old_sign_d  = 1;
ii          = 0;
no_sign_chg = 1;
nsc         = 0;

while (E <= E_max) & no_sign_chg 
    E  = E + E_step;
    ii = ii+1;
    
    % The transfer matrix
    tf_mat       = TransferMatrix(E, k_t);
    c_min(ii)    = -(tf_mat(1,3)*tf_mat(1,4) + tf_mat(2,3)*tf_mat(2,4))/(tf_mat(1,4)^2 + tf_mat(2,4)^2);
    diff_inf(ii) = (tf_mat(1,3) + tf_mat(1,4)*c_min(ii))^2 + (tf_mat(2,3) + tf_mat(2,4)*c_min(ii))^2;
    diff_h(ii)   = tf_mat(3,3) + c_min(ii)*tf_mat(3,4);
    diff_l(ii)   = tf_mat(4,3) + c_min(ii)*tf_mat(4,4);
    
    if (ii > 1)
        
       diff(ii)    = diff_inf(ii) - diff_inf(ii-1);
       new_sign_d  = sign(diff(ii));
       
       if (ii>2)
          if (new_sign_d-old_sign_d==2)
            nsc = nsc+1;
            if (nsc==min_num)
               no_sign_chg = 0; 
            end
          end
       end
       
       old_sign_d = new_sign_d;
    end
end   % -- end of while loop

% Checking if a minimum was found
if (no_sign_chg)    % <-- no minimum
    a     = c_min(ii);
    E     = E;
    error = 1;    
else                % <-- minimum reached
    if (E_step>E_acc)
        % Recursive call, refining the search grid with every step 
        [a,E,error] =  WavefunctionMinimum(E-2*E_step,E,k_t,E_acc); 
    else
        a = c_min(ii-1);
        E = E - E_step;
        error = 0;
    end
end