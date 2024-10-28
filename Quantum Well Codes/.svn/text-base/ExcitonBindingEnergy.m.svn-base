function [beta,EX0,EX0_l] = ExcitonBindingEnergy(z_grid, wf_e, wf_h, vb_index, cb_index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function calculates the exciton binding energy from user
%  supplied wavefunctions (correlated or uncorrelated).
%
%  -> Based on chapter 6 of Harrison. <-
%
%  Theory:
%  - Exciton energy: E = (A+B+C)/D
%  - Exciton binding energy: E = E_e+E_h+E_X0;
%  - The e-h interaction wavefunction: psi_r = exp(-r'/delta)
%  - r' = (x_e-x_h)^2 + (y_e-y_h)^2 + zeta^2*(z_e-z_h)^2
%    zeta=0 - 2D exciton, zeta=1 - 3D exciton
%
%  Input:
%         - wf_e    - electron wavefunction versus z of subband X
%         - wf_h    - hole wavefunction versus z of subband X
%
%  Output:
%         - A       - electron (1-particle) hamiltonian expectation value (<psi|H_e|psi>)                 
%         - B       - hole (1-particle) hamiltonian expectation value (<psi|H_h|psi>)
%         - C       - e-h interaction hamiltonian expectation value (<psi|H_e-h|psi>)
%         - beta    - minimum of Eb (binding energy) vs. lambda (Bohr's radius)
%         - EX0     - minimum binding energy, corresponding lambda
%                     and beta
%         - EX0_l   - minimum binding energy vs. lambda
%         - p       - uncorrelated probability of e-h separation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;
Constants;

%% 1. Definitions

% Dimensionality parameter  
beta_start = 0.001;
beta_step = 0.05;
beta_stop = -1.0;

Consts.epsilon = 13.18*Consts.eps_0;   % relative permittivity of material

% Bohr radius   
lambda_start = 70e-10;       % m
lambda_step = 1e-10;         % m
lambda_stop = -1e-10;        % m

Consts.m = [0.067*Consts.m_0, 0.62*Consts.m_0];
wf.wf_e = wf_e(cb_index,:);
wf.wf_h = wf_h(vb_index,:);

repeat_flag_beta = true;     % repeat variational beta loop flag
repeat_flag_lambda = true;   % repeat variational lambda loop

N_x = 100;

%% 2. Calculations

n = length(wf_e);
d_z = (z_grid(2)-z_grid(1));      % z separation of input potentials
d_a = d_z;                              % separation of adjacent a values 
Eb_min = 1*Consts.e_0;                  % minimum Eb for lambda variation (1eV)
m_xy = Consts.m;                        % e and h x-y plane masses - we assume isotropic mass 
mu_xy = 1/(1/m_xy(1) + 1/m_xy(2));      % exciton reduced mass in x-y plane
pP = CalculateProbabilities(d_a,n,wf);
 
lambda = lambda_start;  % initial Bohr radius 
EX0_l = []; beta_mat = [];

while(1)
    beta = beta_start;
    Eb_min_beta = 1*Consts.e_0;  % minimum Eb for beta variation (1eV)    
   
    while(1)
        [Eb,ABC] = Eb_1S(wf,pP,beta,d_a,lambda,mu_xy,N_x,n);
        
        if (Eb<Eb_min_beta)
            Eb_min_beta = Eb;
            beta_0_lambda = beta;
            repeat_flag_beta = true;
        else
            repeat_flag_beta = false;
        end

        beta = beta + beta_step;
        
        if ((repeat_flag_beta && (beta_stop<0)) || (beta<beta_stop))
            break;
        end
    end

    EX0_l = [EX0_l; lambda/1e-10, Eb_min_beta/(1e-3*Consts.e_0)];
    beta_mat = [beta_mat; lambda/1e-10, beta_0_lambda];
    
    if (Eb_min_beta<Eb_min)
        Eb_min = Eb_min_beta;
        lambda_0 = lambda;
        beta_0 = beta_0_lambda;
        repeat_flag_lambda = true;
    else
        repeat_flag_lambda = false;
    end
    
    lambda = lambda + lambda_step;   % increment Bohr radius

    if ((repeat_flag_lambda && (lambda_stop<0)) || (lambda<lambda_stop))
        Eb_min
        break;
    end
end

EX0 = [Eb_min/(1e-3*Consts.e_0), lambda_0/1e-10, beta_0];

%--------------------------------------------------------------
%---------------- Inner Functions -----------------------------
%--------------------------------------------------------------

function pP = CalculateProbabilities(d_a,n,wf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function calculates the probabilities known as p(a), Pm(a)
%  and Pmu(a) and returns the appropriate structure
%  Uncorrelated probability of finding the electron and hole
%  separated by a distance a.
%
%  Input:
%         - delta_a     - distance between adjacent points
%         - n           - grid size
%         - wf          - wavefucntions structure
%
%  Output:
%         - pP          - structure containing the probabilities vs. a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;

temp_sum = 0;
index_start = 1;

for (ii=0:n-1)
    pP.a(ii+1) = d_a*ii;
    pP.p(ii+1) = 0;
    index = index_start;    % reset wavefunction pointer
    for (jj=0:n-ii-1)
        pP.p(ii+1) = pP.p(ii+1) + ((wf.wf_e(index+ii).^2)*(wf.wf_h(index).^2) + ...
                               (wf.wf_e(index).^2)*(wf.wf_h(index+ii).^2))*d_a;
        index = index + 1;
    end
    temp_sum = temp_sum + pP.p(ii+1)*d_a;
end

function [Eb,ABC] = Eb_1S(wf,pP,beta,d_a,lambda,mu_xy,N_x,n)

global Consts;

A = 0;      % single e hamiltonian expectation value
B = 0;      % single h hamiltonian expectation value
Ct = 0;     % kinetic energy component of C
Cv = 0;     % potential energy component of C
D = 0;      % <psi|psi>
O = 0;      % overlap integral

for (i_a=1:n)
   A = A + pP.p(i_a)*G(pP.a(i_a),beta,lambda,N_x)*d_a;
   B = B + pP.p(i_a)*G(pP.a(i_a),beta,lambda,N_x)*d_a; 
   Ct = Ct + pP.p(i_a)*J(pP.a(i_a),beta,lambda,N_x)*d_a;
   Cv = Cv + pP.p(i_a)*K(pP.a(i_a),beta,lambda,N_x)*d_a;
   D = D + pP.p(i_a)*F(pP.a(i_a),beta,lambda)*d_a;
   
   O = O + wf.wf_e(i_a).*wf.wf_h(i_a)*d_a;
end

A = A*(Consts.hbar^2/(2*Consts.m(1)));
B = B*(Consts.hbar^2/(2*Consts.m(2)));
Ct = Ct*(-Consts.hbar^2/(2*mu_xy));
Cv = Cv*(-Consts.e_0^2/(4*pi*Consts.epsilon));

Eb = (A+B+Ct+Cv)/D;
ABC = [lambda/1e-10, Eb/(1e-3*Consts.e_0), (A+B+Ct)/D/(1e-3*Consts.e_0), Cv/D/(1e-3*Consts.e_0), (O^2)/D]

function g = G(a,beta,lambda,N_x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function returns the value of G(a), to overcome the problem of 
%  divergence when x=0, the integration is performed using a midpoint 
%  sum, the strip width being delta_x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;

delta_x = (1-0)/N_x;
g = 0;

for (x = delta_x/2:delta_x:1)
   g = g + exp(-sqrt(1-(beta^2))*a*(1/x+x)/lambda)*((1-x^2)/(x*(1+x^2)))*delta_x;
end

g = g*(2*pi*((1-(beta^2))^2)*(a^2)/(lambda^2));

function f = F(a,beta,lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function returns the value of F(a).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;

f = 2*pi*lambda*(sqrt(1-(beta^2))*a/2+lambda/4)*exp(-2*sqrt(1-(beta^2))*a/lambda);

function j = J(a,beta,lambda,N_x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function returns the value of J(a), to overcome the problem of 
%  divergence when x=0, the integration is performed using a midpoint 
%  sum, the strip width being delta_x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;

delta_x = (1-0)/N_x;
j13 = 2*pi*(sqrt(1-(beta^2))*a/(2*lambda)-0.25)*exp(-2*sqrt(1-(beta^2))*a/lambda);
j24 = 0;

for (x = delta_x/2:delta_x:1)
   j24 = j24 + (-1/(lambda*((1/x+x)^2)/4) - sqrt(1-(beta^2))*a/((lambda^2)*(1/x+x)/2))*...
               exp(-sqrt(1-(beta^2))*a*(1/x+x)/lambda)*(1/(x^2)-1)*delta_x;        
end

j24 = j24*2*pi*sqrt(1-(beta^2))*a/2;
j = j13 + j24;

function k = K(a,beta,lambda,N_x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function returns the value of K(a), to overcome the problem of 
%  divergence when x=0, the integration is performed using a midpoint 
%  sum, the strip width being delta_x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Consts;

upper_limit = (1-sqrt(1-(beta^2)))/beta;
lower_limit = 0;

delta_x = (upper_limit-lower_limit)/N_x;
k = 0;

for (x = lower_limit+delta_x/2:delta_x:upper_limit)
    k = k + exp(-beta*a*(1/x-x)/lambda)*(1/(x^2)-1)*delta_x;
end

k = k*(2*pi*beta*a/2);
