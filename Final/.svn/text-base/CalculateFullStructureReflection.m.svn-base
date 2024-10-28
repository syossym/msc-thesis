function MCParams = CalculateFullStructureReflection(Structure, Bands, QWParams, Params)

global Consts;

%% Definitions

FullStructureParams = [];
E_vec = min(Params.E_exc):1e-5*Consts.e_0:max(Params.E_exc); % [1.515:0.00005:1.55]*Consts.e_0;
lambda_vec = (2*pi./E_vec)*Consts.c*Consts.hbar; %[1300:0.1:1600].*1e-9;  % [m]
grid = 300;
E_in = [1;0];     % left-hand side incident field amplitude
n_c_l = 1; %M_GaAs.n; % refractive index of the left cladding
n_c_r = 1; %M_GaAs.n; % refractive index of the right cladding

%% Creating structure

delta = Params.delta;

h = waitbar(0,'Building Structure...');
for (ss=1:length(Structure.AllLayers))
    waitbar(ss/length(Structure.AllLayers),h);
    temp_params.x = Structure.AllLayers{ss}.x;
    temp_params.E = 1.525; %(E_vec(1)+0.5*(E_vec(end)-E_vec(1)))/Consts.e_0;
    temp_params.T = Params.T;
    mat = GetMaterial(Structure.AllLayers{ss}.Name, temp_params);
    n_vec_profile(ss) = mat.n;
    n_vec_calc_init(ss,:) = CalculateRefractiveIndex('exp_1', Structure.AllLayers{ss}.x, temp_params.T, E_vec./Consts.e_0);
                            %GetRefractiveIndex(Structure.AllLayers{ss}.Name, E_vec./Consts.e_0);
    l_vec(ss) = delta*Structure.AllLayers{ss}.L*1e-10;      % [m]
end
close(h);

z_vec = 0;
z_grid_calc = 0;
n_profile = n_vec_profile(1);
for (ll=1:length(l_vec))
    z_vec = [z_vec, z_vec(end)+l_vec(ll)];
    vec = linspace(z_grid_calc(end), z_grid_calc(end)+l_vec(ll), grid);
    z_grid_calc = [z_grid_calc, vec ];
    n_profile = [n_profile, ones(1,length(vec))*n_vec_profile(ll)];
end
z_grid_calc = z_grid_calc(2:end);
n_profile = n_profile(2:end);

MCParams.z_grid = z_grid_calc;
MCParams.n_profile = n_profile;

%% Reflection spectrum - Bare MC

h = waitbar(0,'Calculating MC Spectrum');
for (gg=1:length(lambda_vec))
    waitbar(gg/length(lambda_vec),h);

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
close(h);

% figure(11);
% surf(z_vec/1e-10, E_vec/Consts.e_0, I);
% ylabel('E (eV)'); xlabel('z (Angstrom)'); zlabel('|E^{+}+E^{-}|^2');
% title(['\delta=' num2str(delta)]);
% drawnow;

MCParams.BareMC.delta = delta;
MCParams.BareMC.r_MC = r_MC;

%% Reflection spectrum - MC+QW - detuning - HF model

num_mins = Params.num_mins;
N_DEG = Params.N_DEG;
delta_vec = Params.delta_vec; %[0.945:0.0002:0.958];

temp_params.x = 1;
temp_params.E = E_vec/Consts.e_0;
temp_params.T = Params.T;
QW_Mat = GetMaterial('GaAs', temp_params);
Xi_QW = (QW_Mat.n.^2-1)/(4*pi);
Xi_e_Re_TE = interp1(Params.E_exc, (real(QWParams.Xi.TE.HF)), E_vec, 'pchip');
Xi_e_Im_TE = interp1(Params.E_exc, (imag(QWParams.Xi.TE.HF)), E_vec, 'pchip');
Xi_e_Re_TM = interp1(Params.E_exc, (real(QWParams.Xi.TM.HF)), E_vec, 'pchip');
Xi_e_Im_TM = interp1(Params.E_exc, (imag(QWParams.Xi.TM.HF)), E_vec, 'pchip');
Xi_QW_TE = Xi_QW + 1i*Xi_e_Im_TE + Xi_e_Re_TE; 
Xi_QW_TM = Xi_QW + 1i*Xi_e_Im_TM + Xi_e_Re_TM;
n_vec_calc_TE = n_vec_calc_init;
n_vec_calc_TM = n_vec_calc_init;

for (ss=1:length(Structure.AllLayers))
    if (Structure.AllLayers{ss}.IsActive)
        n_QW_TE = smooth(sqrt(1+4*pi*Xi_QW_TE),20);
        n_QW_TM = smooth(sqrt(1+4*pi*Xi_QW_TM),20);
        n_vec_calc_TE(ss,:) = smooth(sqrt(1+4*pi*Xi_QW_TE),20);
        n_vec_calc_TM(ss,:) = smooth(sqrt(1+4*pi*Xi_QW_TM),20);
    end
end

l_vec_init = l_vec/delta;
E_min_r_MCQW = zeros(num_mins, length(delta_vec));
for (dd=1:length(delta_vec))
    %delta_vec(dd)
    l_vec = l_vec_init.*delta_vec(dd);

    %h = waitbar(0,['Calculating MC Spectrum, \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_init(:,gg)];
        %n_vec = [n_c_l, n_vec_profile];
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                                   (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                                (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
                            
        E = [-M_DBR(2,1)/M_DBR(2,2) ; det(M_DBR)/M_DBR(2,2)];

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
        r_MC_detuning(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);
    
    [min_r_MC, min_r_MC_i] = min(abs(r_MC_detuning(dd,:)).^2);
    E_min_r_MC(dd) = E_vec(min_r_MC_i)/Consts.e_0;

    %h = waitbar(0,['Calculating MCQW Spectrum (TE), \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_TE(:,gg)];
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
        r_MCQW_detuning_TE(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);

    %h = waitbar(0,['Calculating MCQW Spectrum (TM), \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_TM(:,gg)];
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
        r_MCQW_detuning_TM(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);

    [maxs_r_MCQW_TE,mins_r_MCQW_TE] = peakdet(smooth(abs(r_MCQW_detuning_TE(dd,:)).^2), 0.001, E_vec/Consts.e_0);
    if (length(mins_r_MCQW_TE) > 0)
        E_min_r_MCQW_TE(1:length(mins_r_MCQW_TE(:,1)),dd) = mins_r_MCQW_TE(:,1);
    else
        E_min_r_MCQW_TE(:,dd) = 0;
    end
    [maxs_r_MCQW_TM,mins_r_MCQW_TM] = peakdet(smooth(abs(r_MCQW_detuning_TM(dd,:)).^2), 0.001, E_vec/Consts.e_0);
    if (length(mins_r_MCQW_TM) > 0)
        E_min_r_MCQW_TM(1:length(mins_r_MCQW_TM(:,1)),dd) = mins_r_MCQW_TM(:,1);
    else
        E_min_r_MCQW_TM(:,dd) = 0;
    end
end

MCParams.N_DEG = N_DEG;
MCParams.Detuning.HF.n_QW_TE = n_QW_TE;
MCParams.Detuning.HF.n_QW_TM = n_QW_TM;
MCParams.Detuning.HF.delta_vec = delta_vec;
MCParams.Detuning.HF.E_vec = E_vec;
MCParams.Detuning.HF.r_MC = r_MC_detuning;
MCParams.Detuning.HF.E_min_MC = E_min_r_MC;
MCParams.Detuning.HF.r_MCQW_TE = r_MCQW_detuning_TE;
MCParams.Detuning.HF.r_MCQW_TM = r_MCQW_detuning_TM;
MCParams.Detuning.HF.E_min_MCQW_TE = E_min_r_MCQW_TE;
MCParams.Detuning.HF.E_min_MCQW_TM = E_min_r_MCQW_TM;

% Plotting
% figure(5);
% for (dd=1:1:length(delta_vec))
%     [maxs_r,mins_r] = peakdet(abs(r_MCQW_detuning(dd,:)).^2, 0.001, E_vec/Consts.e_0);
%     subplot(length(delta_vec),1,dd);
%     plot(E_vec/Consts.e_0, abs(r_MCQW_detuning(dd,:)).^2, 'b', E_vec/Consts.e_0, abs(r_MC_detuning(dd,:)).^2, 'r');
%     axis([min(E_vec)/Consts.e_0, max(E_vec)/Consts.e_0, 0.25, 1]);
%     ylabel(num2str(delta_vec(dd)), 'FontSize', 7);
%     for (mn=1:length(mins_r(:,1)))
%         text(mins_r(mn,1), mins_r(mn,2), [num2str(mins_r(mn,1)) 'eV']);
%     end
%     h_curr = gca;
%     if (dd==1)
%         title_text = ['|r_{MC}|^2 (red), |r_{MCQW}|^2 (blue): T=' num2str(Params.T) 'K, N_{DEG}=' strrep(num2str(N_DEG, '%1.0e'), 'e+0', 'x10^{') '}cm^{-2}'];
%         title(title_text);
%     end
%     if (dd~=length(delta_vec))
%         set(h_curr, 'XTickLabel', '');
%     end
% end
% xlabel('E [eV]');
% 
% figure(6); box on; hold on;
% for (dd=1:1:length(delta_vec))
%     plot(E_vec/Consts.e_0, abs(r_MCQW_detuning(dd,:)).^2, 'b');
%     [min_E_MCQW, min_E_MCQW_i] = min(abs(r_MCQW_detuning(dd,:)).^2);
%     text(E_vec(min_E_MCQW_i)/Consts.e_0, min_E_MCQW, num2str(delta_vec(dd)), 'FontSize', 7);
% end
% ylabel('|r_{MCQW}|.^2'); xlabel('E [eV]');

%% Reflection spectrum - MC+QW - detuning - FCT model

Xi_QW = (QW_Mat.n.^2-1)/(4*pi);
Xi_e_Re_TE = interp1(Params.E_exc, (real(QWParams.Xi.TE.FCT)), E_vec, 'pchip');
Xi_e_Im_TE = interp1(Params.E_exc, (imag(QWParams.Xi.TE.FCT)), E_vec, 'pchip');
Xi_e_Re_TM = interp1(Params.E_exc, (real(QWParams.Xi.TM.FCT)), E_vec, 'pchip');
Xi_e_Im_TM = interp1(Params.E_exc, (imag(QWParams.Xi.TM.FCT)), E_vec, 'pchip');
Xi_QW_TE = Xi_QW + 1i*Xi_e_Im_TE + Xi_e_Re_TE; 
Xi_QW_TM = Xi_QW + 1i*Xi_e_Im_TM + Xi_e_Re_TM;
n_vec_calc_TE = n_vec_calc_init;
n_vec_calc_TM = n_vec_calc_init;

for (ss=1:length(Structure.AllLayers))
    if (Structure.AllLayers{ss}.IsActive)
        n_QW_TE = smooth(sqrt(1+4*pi*Xi_QW_TE),20);
        n_QW_TM = smooth(sqrt(1+4*pi*Xi_QW_TM),20);
        n_vec_calc_TE(ss,:) = smooth(sqrt(1+4*pi*Xi_QW_TE),20);
        n_vec_calc_TM(ss,:) = smooth(sqrt(1+4*pi*Xi_QW_TM),20);
    end
end

l_vec_init = l_vec/delta;
E_min_r_MCQW = zeros(num_mins, length(delta_vec));
for (dd=1:length(delta_vec))
    %delta_vec(dd)
    l_vec = l_vec_init.*delta_vec(dd);

    %h = waitbar(0,['Calculating MC Spectrum, \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_init(:,gg)];
        %n_vec = [n_c_l, n_vec_profile];
        M_DBR = eye(2);
        M_c_l = (1/(2*n_vec(1))).*[(n_vec(1)+n_c_l), (n_vec(1)-n_c_l);...
                                   (n_vec(1)-n_c_l), (n_vec(1)+n_c_l)];
        M_c_r = (1/(2*n_c_r)).*[(n_c_r+n_vec(end)), (n_c_r-n_vec(end));...
                                (n_c_r-n_vec(end)), (n_c_r+n_vec(end))];
                            
        E = [-M_DBR(2,1)/M_DBR(2,2) ; det(M_DBR)/M_DBR(2,2)];

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
        r_MC_detuning(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);
    
    [min_r_MC, min_r_MC_i] = min(abs(r_MC_detuning(dd,:)).^2);
    E_min_r_MC(dd) = E_vec(min_r_MC_i)/Consts.e_0;

    %h = waitbar(0,['Calculating MCQW Spectrum (TE), \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_TE(:,gg)];
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
        r_MCQW_detuning_TE(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);

    %h = waitbar(0,['Calculating MCQW Spectrum (TM), \delta=' num2str(delta_vec(dd))]);
    for (gg=1:length(lambda_vec))
        %waitbar(gg/length(lambda_vec),h);

        n_vec = [n_c_l; n_vec_calc_TM(:,gg)];
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
        r_MCQW_detuning_TM(dd,gg) = -M_DBR(2,1)/M_DBR(2,2);
    end
    %close(h);

    [maxs_r_MCQW_TE,mins_r_MCQW_TE] = peakdet(smooth(abs(r_MCQW_detuning_TE(dd,:)).^2), 0.001, E_vec/Consts.e_0);
    if (length(mins_r_MCQW_TE) > 0)
        E_min_r_MCQW_TE(1:length(mins_r_MCQW_TE(:,1)),dd) = mins_r_MCQW_TE(:,1);
    else
        E_min_r_MCQW_TE(:,dd) = 0;
    end
    [maxs_r_MCQW_TM,mins_r_MCQW_TM] = peakdet(smooth(abs(r_MCQW_detuning_TM(dd,:)).^2), 0.001, E_vec/Consts.e_0);
    if (length(mins_r_MCQW_TM) > 0)
        E_min_r_MCQW_TM(1:length(mins_r_MCQW_TM(:,1)),dd) = mins_r_MCQW_TM(:,1);
    else
        E_min_r_MCQW_TM(:,dd) = 0;
    end
end

MCParams.Detuning.FCT.n_QW_TE = n_QW_TE;
MCParams.Detuning.FCT.n_QW_TM = n_QW_TM;
MCParams.Detuning.FCT.delta_vec = delta_vec;
MCParams.Detuning.FCT.E_vec = E_vec;
MCParams.Detuning.FCT.r_MC = r_MC_detuning;
MCParams.Detuning.FCT.E_min_MC = E_min_r_MC;
MCParams.Detuning.FCT.r_MCQW_TE = r_MCQW_detuning_TE;
MCParams.Detuning.FCT.r_MCQW_TM = r_MCQW_detuning_TM;
MCParams.Detuning.FCT.E_min_MCQW_TE = E_min_r_MCQW_TE;
MCParams.Detuning.FCT.E_min_MCQW_TM = E_min_r_MCQW_TM;
