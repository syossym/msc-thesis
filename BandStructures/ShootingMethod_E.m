function [E,V_e,V_h,m_eff,con_profile] = ShootingMethod_E(p_type, structure, E_ref, num_state, Np, a)

%
% This fucntion uses the shooting technique to calculate the uncorrelated one
% particle energies of any user supplied potential in GaAs/AlGaAs heterostructures.
%
% Input:
%        - 'p_type'    - 'e' for electron / 'hh' for a heavy hole / 'lh' for a light hole.
%        - 'structure' - a matrix describing the structure of the
%                        heterostructure of the form (for a material Ga(1-x)Al(x)As):
%                          [ layer_thickness(A)   x ]
%        - 'E_ref'     - the reference energy for the calculation.
%        - 'num_state' - number of states to calculate.
%        - 'Np'        - number of monolayers in the heterostructure
%
% Output:
%        - 'E' (eV)      - vector of the state energies.
%        - 'V' (J)       - potential profile for e/hh/lh as a function of the growth axis (z).
%        - 'm_eff' (kg)  - effective mass profile for e/hh/lh as a function of the
%                          growth axis (z).
%        - 'con_profile' - concentration profile for the materials in the
%                        alloys in the structure.
% 
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

global Consts;

%% 1. Parameters calculation

% Concentration profile
z_0      = 0;            % z axis in multiples of 1./N
z_size   = 0;            % actual z axis size (A)
N        = 1e-10/a;      % number of points per Angstrom
index    = 1;

% con_profile=[];
% for (ii=1:length(structure(:,1)))
%     con_profile = [con_profile; linspace(z_0,z_0+structure(ii,1)*1e-10,structure(ii,2))', ones(structure(ii,2),1).*structure(ii,3)];
%     z_0 = z_0 + structure(ii,1)*1e-10;
% end
%
% for (ii=1:length(structure(:,1))-1)
%    con_profile(structure(ii,2),2) = 0.5*(con_profile(structure(ii,2)+1,2)-con_profile(structure(ii,2)-1,2)) + con_profile(structure(ii,2)-1,2);
% end

structure(:,1) = structure(:,1)./1e-10;

for (ii=1:length(structure(:,1)))
    n = 0;
    z_size = z_size + structure(ii,1)*1e-10;
    while ((z_0+n*1e-10/N)*1e14 <= z_size*1e14)
        con_profile(index,:) = [(z_0+n*1e-10/N), structure(ii,3)];
        n = n+1;
        index = index+1;
    end
    z_0 = z_0 + n*1e-10/N;
end

if (index-1<Np)
    for (ii=index:Np)
        con_profile(ii,:) = [(con_profile(ii-1,1)+1e-10/N), con_profile(index-1,2)];
        n = n+1;
        if (ii>=Np)
            break;
        end
    end
end

% for (ii=1:length(structure(:,1)))
%     n = 0;
%     z_size = z_size + structure(ii,2)*a;
%     while ((z_0+n*a)*1e14 < z_size*1e14)
%         con_profile(index,:) = [(z_0+n*a), structure(ii,3)];
%         n = n+1;
%         index = index+1;
%     end
%     z_0 = z_0 + n*a;
% end

% Potential profile for e/hh/lh
% for (ii=1:length(con_profile(:,1)))
%     d_V = 1.247*con_profile(ii,2)*e_0;
%     switch(p_type)
%         case 'e'
%             V(ii,:) = [con_profile(ii,1), 0.67*d_V];
%         case 'hh'
%             V(ii,:) = [con_profile(ii,1), 0.33*d_V];
%         case 'lh'
%             error('Data not defined for Ga(1-x)Al(x)As light-hole');
%     end
% end

for (ii=1:length(con_profile(:,1)))
    d_V = 1.247*con_profile(ii,2)*Consts.e_0;
    V_e(ii,:) = [con_profile(ii,1), 0.67*d_V];
    V_h(ii,:) = [con_profile(ii,1), 0.33*d_V];
end

switch(p_type)
    case 'e'
        V = V_e;
    case 'hh'
        V = V_h;
    case 'lh'
        error('Data not defined for Ga(1-x)Al(x)As light-hole');
end

% Effective mass profile for e/hh/lh
for (ii=1:length(con_profile(:,1)))
    switch(p_type)
        case 'e'
            m_eff(ii,:) = [con_profile(ii,1), (0.067+0.083*con_profile(ii,2))*Consts.m_0];
        case 'hh'
            m_eff(ii,:) = [con_profile(ii,1), (0.62+0.14*con_profile(ii,2))*Consts.m_0];
        case 'lh'
            error('Data not defined for Ga(1-x)Al(x)As light-hole');
    end
end

%% Calculating the energies

delta_E = 1e-3*Consts.e_0;     % small energy separation
d_E     = 1e-5*Consts.e_0;     % inf. energy separation
d_z     = V(2,1)-V(1,1);

if (abs(E_ref)<1e-3*Consts.e_0)
    x = min(V(:,2));
else
    x = E_ref;
end

for (ss=1:num_state)

    % Increment energy-search for f(x)=0
    y_temp2 = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));

    while (1)
        y_temp1 = y_temp2;
        x = x + delta_E;
        y_temp2 = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));
        %disp(['State [' num2str(ss) ']- [y1,y2]=[' num2str(y_temp1) ',' num2str(y_temp2) ']: incrementing to ' num2str(x/e_0) 'eV']);
        if (y_temp1*y_temp2<=0)
            %disp(['State [' num2str(ss) ']- Energy search: energy found']);
            break;
        end
    end

    % Improve estimate using midpoint rule
    x = x - abs(y_temp2)/(abs(y_temp1)+abs(y_temp2))*delta_E;
    %disp(['State [' num2str(ss) ']- Energy search: energy improved to ' num2str(x/Consts.e_0) 'eV']);

    % Newton-Raphson method
    while(1)
        y  = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));
        dy = (WfAtInf(x+d_E,d_z,V,m_eff,length(V(:,1))) - ...
              WfAtInf(x-d_E,d_z,V,m_eff,length(V(:,1)))) / (2*d_E);

        if (dy == 0)
            break;
        end

        x  = x - y/dy;

        if (abs(y/dy)<=1e-12*Consts.e_0)
            break;
        end
    end

    E(ss) = x/Consts.e_0;
    x     = x + delta_E;
    %disp(['State [' num2str(ss) '] energy is ' num2str(E(ss)) 'eV']);

end

V = V;

% ------------ Inner Functions ----------------------

function wf_inf = WfAtInf(E,d_z,V,m_eff,N)

% This function returns the value of the wavefunction
% at +infinity for a given value of the energy. The solution
% to the energy occurs for wf(+infinity)=0.

global Consts;

% Wavefunction at z-delta_z, z and z+delta_z
wf(1) = 0.0;
wf(2) = 1.0;
%disp(['-- E=' num2str(E)]);
for (ii=2:N-1)   % last potential not used
    %disp(['ii=' num2str(ii) ', V-E=' num2str(V(ii,2)-E) ', m=[' num2str(m_eff(ii-1,2)) ',' num2str(m_eff(ii,2)) ',' num2str(m_eff(ii+1,2)) ']']);
    wf(3) = ( ...
                ( ...
                2*((d_z/Consts.hbar)^2)*(V(ii,2)-E) + ...
                2./(m_eff(ii,2)+m_eff(ii+1,2)) + ...
                2./(m_eff(ii,2)+m_eff(ii-1,2)) ...
                )*wf(2) - (2./(m_eff(ii,2)+m_eff(ii-1,2)))*wf(1) ...
                ) * ((m_eff(ii,2)+m_eff(ii+1,2)))/2;
    wf(1) = wf(2);
    wf(2) = wf(3);
    %disp(['wf=[' num2str(wf) ']']);
end

wf_inf = wf(3);