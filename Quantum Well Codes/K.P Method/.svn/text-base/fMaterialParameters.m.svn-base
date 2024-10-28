
function s_Params = fMaterialParameters(v_MaterialName, v_AlFraction)
%----------------------------------------------------------------------		
%   Name   : fMaterialParameters 
%   Type   : Function
%	
%   Description : 										   
%	This function creates a structure containing the 
%   material parameters of the specified material. 
%   All parameters are taken at 300K.
%   
%   Sources:
%        - http://www.ioffe.rssi.ru/SVA/NSM/Semicond/
%		 - Physics of Optoelectronic Devices, Chuang
%
%   Input  : (string) v_MaterialName  - the material name. 
%			 (int) v_AlFraction - Al fraction in the alloy Al(x)Ga(1-x)As. 
%   Output : (struct) s_Params - structure contaning the 
%                                material properties.
%         
%   Tested : Matlab 7.5.0
%	    By : Yossi Michaeli, May 2009
%	
%----------------------------------------------------------------------

global s_PhysConsts;

if (strcmp(v_MaterialName,'AlGaAs') && nargin<2)
	error('fMaterialParameters:Nargin_Error',...
	      'Number of arguments for [AlGaAs] should be 2');
end

s_Params.Name = v_MaterialName;

switch (v_MaterialName)
	case 'GaAs'
	
		s_Params.a    = 5.65325;   				% [Angstrom]
		s_Params.m_e  = 0.063*s_PhysConsts.m_0; % [kg]
		s_Params.m_hh = 0.51*s_PhysConsts.m_0;  % [kg]
		s_Params.m_lh = 0.082*s_PhysConsts.m_0; % [kg]
		s_Params.E_g  = 1.424;                  % [eV]
		s_Params.E_p  = 25.7;                   % Optical matrix parameter [eV]
		s_Params.Del  = 0.34;                   % [eV]
		s_Params.E_v  = -6.92;                  % [eV]
		s_Params.x    = 0;						% the fraction of [Al] in the alloy
		
		% Luttinger params
		s_Params.g1   = 6.85;
		s_Params.g2   = 2.06; 
		s_Params.g3   = 2.93;
	
	case 'AlAs'
	
		s_Params.a    = 5.6600;   				% [Angstrom]
		s_Params.m_e  = 0.15*s_PhysConsts.m_0;  % [kg]
		s_Params.m_hh = 0.79*s_PhysConsts.m_0;  % [kg]
		s_Params.m_lh = 0.15*s_PhysConsts.m_0;  % [kg]
		s_Params.E_g  = 3.03;                   % [eV]
		s_Params.E_p  = 21.1;                   % Optical matrix parameter [eV]
		s_Params.Del  = 0.28;                   % [eV]
		s_Params.E_v  = -7.49;                  % [eV]
		s_Params.x    = 1;						% the fraction of [Al] in the alloy
		
		% Luttinger params
		s_Params.g1   = 3.45;
		s_Params.g2   = 0.68; 
		s_Params.g3   = 1.29;
	
	case 'AlGaAs'
	
		s_Params.a    = 5.6533+0.00788*v_AlFraction;     % [Angstrom]
		s_Params.m_e  = (0.063+0.083*v_AlFraction)*...
									  s_PhysConsts.m_0;  % [kg]
		s_Params.m_hh = (0.51+0.25*v_AlFraction)*...
									  s_PhysConsts.m_0;  % [kg]
		s_Params.m_lh = (0.082+0.068*v_AlFraction)*...
									  s_PhysConsts.m_0;  % [kg]
		s_Params.x    = v_AlFraction;			         % the fraction of [Al] in the alloy
		
		if (v_AlFraction<0.45)
			s_Params.E_g  = 1.424+1.247*v_AlFraction;    % [eV]
		else									  
			s_Params.E_g  = 1.9+0.125*v_AlFraction+...
					        0.143*v_AlFraction^2;        % [eV]
		end
		
		s_Params.E_p  = v_AlFraction*21.1+...
						(1-v_AlFraction)*25.7;           % Optical matrix parameter [eV]
		s_Params.Del  = v_AlFraction*0.68+...
						(1-v_AlFraction)*0.34; 		     % [eV]
		s_Params.E_v  = v_AlFraction*(-7.49)+...
						(1-v_AlFraction)*(-6.92);        % [eV]                 
		
		% Luttinger params
		s_Params.g1   = v_AlFraction*3.45+...
						(1-v_AlFraction)*6.85;
		s_Params.g2   = v_AlFraction*0.68+...
						(1-v_AlFraction)*2.06;
		s_Params.g3   = v_AlFraction*1.29+...
						(1-v_AlFraction)*2.93;
	
	otherwise
		error('fMaterialParameters:Unknown_Material',...
	          ['Unknown material [', num2str(v_MaterialName),']']);
end

v_dV = (1.247*s_Params.x*s_PhysConsts.e_0);
s_Params.V_c = 0.67*v_dV; 
s_Params.V_v = 0.33*v_dV;


% Constract parameters
% --------------------

s_Params.gc1  = (s_PhysConsts.hbar^2)*s_Params.g1/...
				(2*s_PhysConsts.m_0*s_PhysConsts.e_0*(s_Params.a^2));
s_Params.gc2  = (s_PhysConsts.hbar^2)*s_Params.g2/...
				(2*s_PhysConsts.m_0*s_PhysConsts.e_0*(s_Params.a^2));
s_Params.gc3  = (s_PhysConsts.hbar^2)*s_Params.g3/...
				(2*s_PhysConsts.m_0*s_PhysConsts.e_0*(s_Params.a^2));