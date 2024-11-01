function M = GetMaterial(name, params)
%
% This function gets the paterial parameters.
%
%   Input:  'name'   - the name of the material
%           'params' - structure containing parameters of the specified
%                      material (relevant only for ternary materials)
%
%   Output: 'M' - structure containing the following fields:
%                TBD
%
%   Note: For all ternary materials the following interpolation equation
%         is used - P(A(x)B(1-x)C)=xP(AC)+(1-x)P(BC)
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, September 2009
% Edited by: -
%

global Consts;

if ((strcmp(name,'AlGaAs') || strcmp(name,'GaAlAs')) && nargin<2)
    error('GetMaterial:Nargin_Error',...
        'Number of arguments for [AlGaAs] should be 2');
end

M.Name = name;

switch (name)
    
    case 'AlGaAs'  %-------Al(x)Ga(1-x)As--------------------------------
    case 'GaAlAs'
        
        M.a     = (5.6533+0.00788*params.x)*1e-10;     % [m]
        M.m_e   = (0.063+0.083*params.x)*Consts.m_0;   % [kg]
        M.m_hh  = (0.51+0.25*params.x)*Consts.m_0;     % [kg]
        M.m_lh  = (0.082+0.068*params.x)*Consts.m_0;   % [kg]
        M.eps_r = 12.90-2.84*params.x;
        M.x     = params.x;	        		           % the fraction of [Al] in the alloy
        
        if (params.x<0.45)
            M.E_g_300  = 1.424+1.247*params.x;                 % Gamma [eV]
        else
            M.E_g_300  = 1.9+0.125*params.x+0.143*params.x^2;  % X [eV]
            %M.E_g_300_G  = 1.424+1.247*params.x+1.247*params.x^2;  % Gamma [eV]
        end
        %M.E_g_300 = 1.425+1.155*params.x+0.37*params.x^2;
        %M.E_g_300 = (1-params.x)*1.424 + params.x*2.12;
        M.E_g = GetMaterialBandGap(params.x,params.T);
        
        M.E_p  = params.x*21.1+(1-params.x)*25.7;          % Optical matrix parameter [eV]
        M.Del  = params.x*0.68+(1-params.x)*0.34; 	       % [eV]
        M.E_v  = params.x*(-7.49)+(1-params.x)*(-6.92);    % [eV]
        
        if (params.E ~= 0)
            %M.n = params.x*GetRefractiveIndex('AlAs',params.E)...
            %    +(1-params.x)*GetRefractiveIndex('GaAs',params.E);
            M.n = GetRefractiveIndex('GaAlAs',params.E);
        else
            M.n =  params.x*(3)+(1-params.x)*(3.3);
        end
        
        % Luttinger params
        M.g1   = params.x*3.45+(1-params.x)*6.85;
        M.g2   = params.x*0.68+(1-params.x)*2.06;
        M.g3   = params.x*1.29+(1-params.x)*2.93;
        
    case 'GaAs'    %-----------------------------------------------------
        
        M.a    = 5.65325*1e-10;          % [m]
        M.m_e  = 0.063*Consts.m_0;       % Gamma-valley mass [kg]
        M.m_hh = 0.51*Consts.m_0;        % [kg]
        M.m_lh = 0.082*Consts.m_0;       % [kg]
        M.E_g_300 = 1.424;               % at 300K [eV]
        M.E_g_0 = 1.521;                 % at 0K [eV]
        M.E_g = GetMaterialBandGap(1,params.T);
        M.E_p  = 25.7;                   % optical matrix parameter [eV]
        M.Del  = 0.34;                   % [eV]
        M.E_v  = -6.92;                  % [eV]
        M.eps_r = 12.4;
        M.x    = 0;	      	     	     % the fraction of [Al] in the alloy
        if (params.E ~= 0)
            M.n = GetRefractiveIndex('GaAs',params.E);
        else
            M.n = 3.3;
        end
        
        % Luttinger params
        M.g1   = 6.85;
        M.g2   = 2.06;
        M.g3   = 2.93;
        
    case 'AlAs'    %-----------------------------------------------------
        
        M.a    = 5.6600*1e-10;   	     % [m]
        M.m_e  = 0.15*Consts.m_0;        % [kg]
        M.m_hh = 0.79*Consts.m_0;        % [kg]
        M.m_lh = 0.15*Consts.m_0;        % [kg]
        M.E_g_300 = 2.163;               % [eV]
        M.E_g_0 = 2.239;                 % [eV]
        M.E_g = GetMaterialBandGap(1,params.T);  % [eV]
        M.E_p  = 21.1;                   % Optical matrix parameter [eV]
        M.Del  = 0.28;                   % [eV]
        M.E_v  = -7.49;                  % [eV]
        M.eps_r = 10.1;
        M.x    = 1;		    		     % the fraction of [Al] in the alloy
        if (params.E ~= 0)
            M.n = GetRefractiveIndex('AlAs',params.E);
        else
            M.n = 3;
        end
        
        % Luttinger params
        M.g1   = 3.45;
        M.g2   = 0.68;
        M.g3   = 1.29;
        
    otherwise % ---------------------------------------------------------
        error('fMaterialParameters:Unknown_Material',...
            ['Unknown material [', num2str(v_MaterialName),']']);
end


