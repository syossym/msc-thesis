
function [s_StructParams] = fQWStructure(v_AlFraction, v_WellWidth, v_MonolayerSize, v_CladdingWellRatio)
%-------------------------------------------------------------------------------------		
%   Name   : fQWStructure 
%   Type   : Function
%	
%   Description : 										   
%	This function creates the Quantum Well profile parameters.
%   
%
%   Input  : (double) v_WellWidth - the width of the QW in [Angstrom].
%			 (double) v_MonolayerSize - the size of a single monolayer [Angstrom].
%            (int) v_CladdingWellRatio - cladding width / well width.
%		     (struct) s_Materials - structure containing the parameters of the materials.
%   Output : (struct) s_Params - the structure parameters.
%         
%   Tested : Matlab 7.5.0
%	    By : Yossi Michaeli, May 2009
%	
%-------------------------------------------------------------------------------------

global s_PhysConsts;

if (nargin<2)
	error('fQWStructure:Nargin_Error',...
	      'Number of arguments should be 2');
end

% Structure Parameters (Well and Cladding layers)
% -----------------------------------------------

s_StructParams.Well = fMaterialParameters('GaAs');
s_StructParams.Clad = fMaterialParameters('AlGaAs', v_AlFraction);

s_StructParams.Well.W      = v_WellWidth; 				                     % well width [Angstrom]
s_StructParams.Clad.W 	   = v_CladdingWellRatio*s_Params.W_w;               % cladding bulk width [Angstrom]
s_StructParams.Well.N      = round(v_WellWidth/v_MonolayerSize);             % number of well layers (each layer is v_MonolayerSize width)
s_StructParams.Clad.N      = v_CladdingWellRatio*s_Params.N_w;               % number of bulk cladding layers
s_StructParams.N_tot       = 2*s_Materials.Clad.N+s_Materials.Well.N;        % total number of layers

if (s_StructParams.Clad.x < 1 && s_StructParams.Clad.x > 0)
	s_StructParams.Struct  = [s_StructParams.Clad.W, s_StructParams.Clad.N, s_StructParams.Clad.x ; 
					          s_StructParams.Well.W, s_StructParams.Well.N, s_StructParams.Well.x ; 
					          s_StructParams.Clad.W, s_StructParams.Clad.N, s_StructParams.Clad.x];
else
	error('fQWStructure:Bulk_Material_Error',...
	      'The [Al] fraction in the cladding layer material alloy should be in the [0,1] range');
end

% Parameter profiles
% ------------------

v_z_0      = 0;      % z axis
V_z_size   = 0;      % actual z axis size [m]
index    = 1;

for (ii=1:length(s_StructParams.Struct(:,1)))
    n = 0;
    V_z_size = V_z_size + s_StructParams.Struct(ii,1)*1e-10;
    while ((v_z_0+n*v_MonolayerSize*1e-10)*1e14 <= V_z_size*1e14)
        v_Profile(index,:) = [(z_0+n*v_MonolayerSize*1e-10), ii];
        n = n+1;
        index = index+1;
    end
    v_z_0 = v_z_0 + n*v_MonolayerSize*1e-10;
end

if (index-1<s_StructParams.N_tot)
    for (ii=index:s_StructParams.N_tot)
        v_Profile(ii,:) = [(v_Profile(ii-1,1)+v_MonolayerSize*1e-10), v_Profile(index-1,2)];
        n = n+1;
        if (ii>=s_StructParams.N_tot)
            break;
        end
    end
end

s_StructParams.Profile.z = v_Profile(:,1);      % z-axis [m]

v_CladdingIndices = v_Profile(:,2)==1 || v_Profile(:,2)==3;
v_WellIndices = v_Profile(:,2)==2;

s_StructParams.Profile.x                    = v_Profile(:,2);      		    % Al concentration profile
s_StructParams.Profile.x(v_CladdingIndices) = s_StructParams.Clad.x;
s_StructParams.Profile.x(v_WellIndices)     = s_StructParams.Well.x;

s_StructParams.Profile.V_c                    = v_Profile(:,2);    			% conductance band profile [eV]
s_StructParams.Profile.V_c(v_CladdingIndices) = s_StructParams.Clad.V_c;
s_StructParams.Profile.V_c(v_WellIndices)     = s_StructParams.Well.V_c;

s_StructParams.Profile.V_v                    = v_Profile(:,2);    			% valence band profile [eV]
s_StructParams.Profile.V_v(v_CladdingIndices) = s_StructParams.Clad.V_v;
s_StructParams.Profile.V_v(v_WellIndices)     = s_StructParams.Well.V_v;

s_StructParams.Profile.m_e                    = v_Profile(:,2);    			% electron effective mass [kg]
s_StructParams.Profile.m_e(v_CladdingIndices) = s_StructParams.Clad.m_e;
s_StructParams.Profile.m_e(v_WellIndices)     = s_StructParams.Well.m_e;

s_StructParams.Profile.m_hh                    = v_Profile(:,2);    		% heavy hole effective mass [kg]
s_StructParams.Profile.m_hh(v_CladdingIndices) = s_StructParams.Clad.m_hh;
s_StructParams.Profile.m_hh(v_WellIndices)     = s_StructParams.Well.m_hh;

s_StructParams.Profile.m_lh 				   = v_Profile(:,2);    	    % light hole effective mass [kg]
s_StructParams.Profile.m_lh(v_CladdingIndices) = s_StructParams.Clad.m_lh;
s_StructParams.Profile.m_lh(v_WellIndices)     = s_StructParams.Well.m_lh;