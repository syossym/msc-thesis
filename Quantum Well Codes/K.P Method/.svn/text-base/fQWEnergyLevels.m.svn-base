
function [s_EnergyLevels] = fQWEnergyLevels(s_StructParams, v_ParticleType)
%-------------------------------------------------------------------------------------		
%   Name   : fQWEnergyLevels
%   Type   : Function
%	
%   Description : 										   
%	This function calculates the energy levels in the quantum well.
%   
%
%   Input  : 
%   Output : 
%         
%   Tested : Matlab 7.5.0
%	    By : Yossi Michaeli, May 2009
%	
%-------------------------------------------------------------------------------------

global s_PhysConsts;

switch (v_ParticleType)
	case 'e'
		v_m_eff = s_StructParams.Profile.m_e;
		v_V = s_StructParams.Profile.V_c;
	case 'hh'
		v_m_eff = s_StructParams.Profile.m_hh;
		v_V = s_StructParams.Profile.V_v;
	case 'lh'
		v_m_eff = s_StructParams.Profile.m_lh;
		v_V = s_StructParams.Profile.V_v;
end

v_delta_E = 1e-3*s_PhysConsts.e_0;     % small energy separation
v_dE      = 1e-5*s_PhysConsts.e_0;     % inf. energy separation
v_dz      = v_V(2,1)-v_V(1,1);

x = min(V(:,2));

for (ss=1:num_state)

    % Increment energy-search for f(x)=0
    y_temp2 = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));

    while (1)
        y_temp1 = y_temp2;
        x = x + delta_E;
        y_temp2 = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));
        %disp(['State [' num2str(ss) ']- [y1,y2]=[' num2str(y_temp1) ',' num2str(y_temp2) ']: incrementing to ' num2str(x/e_0) 'eV']);
        if (y_temp1*y_temp2<=0)
            disp(['State [' num2str(ss) ']- Energy search: energy found']);
            break;
        end
    end

    % Improve estimate using midpoint rule
    x = x - abs(y_temp2)/(abs(y_temp1)+abs(y_temp2))*delta_E;
    disp(['State [' num2str(ss) ']- Energy search: energy improved to ' num2str(x/e_0) 'eV']);

    % Newton-Raphson method
    while(1)
        y  = WfAtInf(x,d_z,V,m_eff,length(V(:,1)));
        dy = (WfAtInf(x+d_E,d_z,V,m_eff,length(V(:,1))) - ...
            WfAtInf(x-d_E,d_z,V,m_eff,length(V(:,1)))) / (2*d_E);

        if (dy == 0)
            break;
        end

        x  = x - y/dy;

        if (abs(y/dy)<=1e-12*e_0)

            break;
        end
    end

    E(ss) = x/(e_0);
    x     = x + delta_E;
    disp(['State [' num2str(ss) '] energy is ' num2str(E(ss)) 'meV']);

end

V = V;

% ---------------- Inner functions --------------------%

function wf_inf = fWaveFunctionAtInf(E,d_z,V,m_eff,N)

% This function returns the value of the wavefunction
% at +infinity for a given value of the energy.  The solution
% to the energy occurs for wf(+infinity)=0.

global hbar e_0 a m_0;
format long;

% Wavefunction at z-delta_z, z and z+delta_z
wf(1) = 0.0;
wf(2) = 1.0;
%disp(['-- E=' num2str(E)]);
for (ii=2:N-1)   % last potential not used
    %disp(['ii=' num2str(ii) ', V-E=' num2str(V(ii,2)-E) ', m=[' num2str(m_eff(ii-1,2)) ',' num2str(m_eff(ii,2)) ',' num2str(m_eff(ii+1,2)) ']']);
    wf(3) = ( ...
        ( ...
        2*((d_z/hbar)^2)*(V(ii,2)-E) + ...
        2./(m_eff(ii,2)+m_eff(ii+1,2)) + ...
        2./(m_eff(ii,2)+m_eff(ii-1,2)) ...
        )*wf(2) - (2./(m_eff(ii,2)+m_eff(ii-1,2)))*wf(1) ...
        ) * ((m_eff(ii,2)+m_eff(ii+1,2)))/2;
    wf(1) = wf(2);
    wf(2) = wf(3);
    %disp(['wf=[' num2str(wf) ']']);
end

wf_inf = wf(3);