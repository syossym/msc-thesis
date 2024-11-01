function tf_total = TransferMatrix(E, k_t)
% This function returns the transfer matrix for the whole structure 
% for E and k_t. 
%
% Inputs:
% - E - subband energy (J)
% - k_t - in-plane hole momentum

% Global variables
global gamma1  gamma2  gamma3  V  N  dz;

% Check the k_t value
if ~exist('k_t')
    k_t = 0;
end

%% Building the transfer matrix

tf_total = eye(4);
tf       = zeros(4,4);
tf(1,3)  = 1; 
tf(2,4)  = 1;

for (n=2:N-1)
  
    if (gamma1(n)~=gamma1(n-1) | n==2)
        
        C1       = gamma1(n)^2 - 4*gamma2(n)^2;
        Cf       = 1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1;
        
        tf(3,1)  = (-1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1) / Cf;
        tf(3,2)  = (2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)-2*gamma2(n))) / Cf;
        
        tf_33    = ( 2 + (k_t^2)*(dz^2)*(gamma1(n)+gamma2(n))/(gamma1(n)-2*gamma2(n)) - 3*(k_t^3)*(dz^3)*gamma2(n)*gamma3(n)/C1 ) / Cf;
        tf_34    = ( sqrt(3)*(k_t^2)*(dz^2)*gamma2(n)/(gamma1(n)-2*gamma2(n)) - 2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)-2*gamma2(n)) - sqrt(3)*(k_t^3)*(dz^3)*gamma3(n)*(gamma1(n)-gamma2(n))/C1 ) / Cf;
   
        tf(4,1)  = (-2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)+2*gamma2(n))) / Cf;
        tf(4,2)  = (-1 + 3*(k_t^2)*(dz^2)*(gamma3(n)^2)/C1) / Cf;
        
        tf_43    = ( sqrt(3)*(k_t^2)*(dz^2)*gamma2(n)/(gamma1(n)+2*gamma2(n)) + 2*sqrt(3)*k_t*dz*gamma3(n)/(gamma1(n)+2*gamma2(n)) + sqrt(3)*(k_t^3)*(dz^3)*gamma3(n)*(gamma1(n)+gamma2(n))/C1) / Cf;
        tf_44    = ( 2 + (k_t^2)*(dz^2)*(gamma1(n)-gamma2(n))/(gamma1(n)+2*gamma2(n)) + 3*(k_t^3)*(dz^3)*gamma2(n)*gamma3(n)/C1) / Cf;
        
        gammap   = gamma1(n) + 2*gamma2(n);
        gammam   = gamma1(n) - 2*gamma2(n);
        const34  = sqrt(3)*k_t*(dz^3)*gamma3(n)/C1;
        
    end
    
    dE       = V(n) - E;
    
    tf(3,3)  = tf_33 + ( (dz^2)*dE/gammam ) / Cf; 
    tf(3,4)  = tf_34 + ( -const34*dE ) / Cf;
    tf(4,3)  = tf_43 + ( const34*dE ) / Cf;
    tf(4,4)  = tf_44 + ( (dz^2)*dE/gammap ) / Cf;
    
    tf_total = tf * tf_total;
    
    if (n==N-10)
       tf_p = tf_total; 
    end
end
