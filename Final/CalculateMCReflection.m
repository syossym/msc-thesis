function MCParams = CalculateMCReflection(Structure, Params, width_factor)

global Consts;

if (nargin < 3)
    width_param = 1;
end

%% Definitions

FullStructureParams = [];
E_vec = Params.E_exc;
lambda_vec = (2*pi./E_vec)*Consts.c*Consts.hbar; %[1300:0.1:1600].*1e-9;  % [m]
grid = 300;
E_in = [1;0];     % left-hand side incident field amplitude
n_c_l = 1; %M_GaAs.n; % refractive index of the left cladding
n_c_r = 1; %M_GaAs.n; % refractive index of the right cladding

%% Creating structure

delta = Params.delta;

%h = waitbar(0,'Building Structure...');
for (ss=1:length(Structure.AllLayers))
    %waitbar(ss/length(Structure.AllLayers),h);
    temp_params.x = Structure.AllLayers{ss}.x;
    temp_params.E = 1.525; %(E_vec(1)+0.5*(E_vec(end)-E_vec(1)))/Consts.e_0;
    temp_params.T = Params.T;
    mat = GetMaterial(Structure.AllLayers{ss}.Name, temp_params);
    n_vec_profile(ss) = mat.n;
    n_vec_calc_init(ss,:) = CalculateRefractiveIndex('exp_1', Structure.AllLayers{ss}.x, temp_params.T, E_vec./Consts.e_0);
    %ModifiedSellmeierEquation(temp_params.x, temp_params.T, E_vec); 
    %GetRefractiveIndex(Structure.AllLayers{ss}.Name, E_vec./Consts.e_0, Structure.AllLayers{ss}.x); % 
    l_vec(ss) = delta*Structure.AllLayers{ss}.L*1e-10*width_factor;      % [m]
end
%close(h);

z_vec = 0;
z_grid_calc = 0;
n_profile = n_vec_calc_init(1,1);
for (ll=1:length(l_vec))
    z_vec = [z_vec, z_vec(end)+l_vec(ll)];
    vec = linspace(z_grid_calc(end), z_grid_calc(end)+l_vec(ll), grid);
    z_grid_calc = [z_grid_calc, vec ];
    n_profile = [n_profile, ones(1,length(vec))*n_vec_calc_init(ll,1)];
end
z_grid_calc = z_grid_calc(2:end);
n_profile = n_profile(2:end);

MCParams.z_grid = z_grid_calc;
MCParams.n_profile = n_profile;

% figure('Name', 'Refractive index profile'); 
% subplot(211); plot(z_grid_calc/1e-10, real(n_profile));
% ylabel('\Re(n)');
% subplot(212); plot(z_grid_calc/1e-10, imag(n_profile));
% ylabel('\Im(n)'); xlabel('z (Angstrom)');

%% Reflection spectrum - Bare MC

%h = waitbar(0,'Calculating MC Spectrum');
for (gg=1:length(lambda_vec))
    %waitbar(gg/length(lambda_vec),h);

    n_vec = [n_c_l; n_vec_calc_init(:,gg)];
    %n_vec = [n_c_l, n_vec_profile];
    M_DBR = eye(2);
    M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
        (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
    M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
        (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
    %M_DBR = M_c_l*M_DBR;
    %phi_vec(gg) = (2*pi/lambda_vec(gg))*(n_vec(1)*l_vec(1));
    for (nn=1:length(l_vec))
        k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
        M_p = [exp(1i*k*l_vec(nn))   ,         0          ;...
            0           , exp(-1i*k*l_vec(nn))];
        M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
            (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];

        M_DBR = M_p*M_i*M_DBR;
    end
    %M_DBR = M_c_r*M_DBR;
    r_MC(gg) = -M_DBR(2,1)/M_DBR(2,2);
    
%     % Field intensity
%     E_in = [1; r_MC(gg)];
%     z_vec(1) = 0; index = 1;
%     I_temp = [];
%     M_DBR = eye(2);
%     for (nn=1:length(l_vec))
%         k = (2*pi/lambda_vec(gg))*n_vec(nn+1);
%         M_i = (1/(2*n_vec(nn+1))).*[(n_vec(nn+1)+n_vec(nn)) , (n_vec(nn+1)-n_vec(nn)) ; ...
%             (n_vec(nn+1)-n_vec(nn)) , (n_vec(nn+1)+n_vec(nn))];
%         dz = l_vec(nn)/10;
%         M_DBR = M_i*M_DBR;
%         
%         for (zz = 1:10)  
%             M_p = [exp(1i*k*dz)   ,  0  ;...
%                          0  , exp(-1i*k*dz)];
%             M_DBR = M_p*M_DBR;
%             E = M_DBR*E_in;
%             I_temp = [I_temp, (E(1)+E(2))*conj(E(1)+E(2))];
%             index = index + 1;
%             z_vec(index) = z_vec(index-1)+dz;
%         end
%     end
%     I(gg,:) = I_temp;
end
%close(h);

% figure(11);
% surf(z_vec/1e-10, E_vec/Consts.e_0, I);
% ylabel('E (eV)'); xlabel('z (Angstrom)'); zlabel('|E^{+}+E^{-}|^2');
% title(['\delta=' num2str(delta)]);
% drawnow;

MCParams.delta = delta;
MCParams.r_MC = r_MC;
MCParams.E_vec = Params.E_exc;
