function wfs = OneParticleWavefunctions(p_type,E,V,m_eff,con_profile)

% This program uses a shooting technique to calculate the
% uncorrelated one particle energies of any user supplied
% potential.

global hbar e_0 a m_0; 

for (ii=1:length(E))
    [wf,N] = CalcWf(E(ii)*e_0,V,m_eff);
    wfs(ii,:) = wf/sqrt(N);
end

function [wf,N] = CalcWf(E,V,m_eff)

global hbar e_0 a m_0; 

d_z = V(2,1) - V(1,1);
N   = 0;                  % normalization integral

% Boundary conditions
psi(1) = 0.0; wf(1) = psi(1);
psi(2) = 1.0; wf(2) = psi(2);

N = N + psi(1)^2; N = N + psi(2)^2;

for (ii=2:length(V(:,1))-1)
    psi(3) = ( ...
              ( ...
                  2*((d_z/hbar)^2)*(V(ii,2)-E) + ...
                  2./(m_eff(ii,2)+m_eff(ii+1,2)) + ...
                  2./(m_eff(ii,2)+m_eff(ii-1,2)) ...
              )*psi(2) - (2./(m_eff(ii,2)+m_eff(ii-1,2)))*psi(1) ...
             ) * (m_eff(ii,2)+m_eff(ii+1,2))/2;

    wf(ii+1) = psi(3);
    N  = N + psi(3)^2;
    psi(1) = psi(2);
    psi(2) = psi(3);
end

N = N*d_z;



