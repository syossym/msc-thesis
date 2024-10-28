function WFs = ShootingMethod_WFs(p_type,E,V,m_eff,con_profile)

%
% This fucntion uses the shooting technique  to calculate the
% uncorrelated one particle energies of any user supplied
% potential.
%
% Input:
%        - 'p_type'    - 'e' for electron / 'hh' for a heavy hole / 'lh'
%                         for a light hole
%        - 'E'         - the level energy
%        - 'V'         - potential profile for e/hh/lh as a function of the
%                        growth axis (z)
%        - 'm_eff'     - effective mass profile for e/hh/lh as a function of the
%                        growth axis (z)
%        - 'con_profile' - concentration profile for the materials in the
%                        alloys in the structure
%
% Output:
%        - 'WFs'   - vector of the wavefunction for the specified level as
%                    a fucntion of the growth axis (z)
%         
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

global Consts; 

for (ii=1:length(E))
    [wf,N] = CalcWf(E(ii)*Consts.e_0,V,m_eff);
    WFs(ii,:) = wf/sqrt(N);
end

% -------------- Inner Functions ------------------------

function [wf,N] = CalcWf(E,V,m_eff)

global Consts; 

d_z = V(2,1) - V(1,1);
N   = 0;                  % normalization integral

% Boundary conditions
psi(1) = 0.0; wf(1) = psi(1);
psi(2) = 1.0; wf(2) = psi(2);

N = N + psi(1)^2; N = N + psi(2)^2;

for (ii=2:length(V(:,1))-1)
    psi(3) = ( ...
              ( ...
                  2*((d_z/Consts.hbar)^2)*(V(ii,2)-E) + ...
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



