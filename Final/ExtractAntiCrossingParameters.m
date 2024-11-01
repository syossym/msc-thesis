function ACParams = ExtractAntiCrossingParameters(MCParams, Fits, Params)

global Consts;

close all; ACParams = [];

[n,m] = size(MCParams);
if (m==1 && n==1)
    dim = '(1,1)';
    index = '{1,1}';
elseif (m==1)
    dim = '(:,1)';
    index = '{ii,1}';
elseif(n==1)
    dim = '(1,:)';
    index = '{1,ii}';
end

% test_vec = 1.52:0.0000001:1.54;     % [eV]

%% Coupled Oscillator Model With the Calculated parameters

try
    disp(' - Calculating the anti-crossing curves using the extracted parameters and the coupled oscillator model');
    for (ii=1:eval(['length(MCParams' dim ' )']))
        eval(['current_vec = MCParams' index ';']);
        
        % Get the model parameters
        width_index = round(length(current_vec.Detuning.delta_vec)/2);
        gamma_MC = GetMCLineWidth(current_vec);
        [gamma_X_TE, gamma_X_TM] = GetExcitonLineWidths(current_vec, width_index);
        
        % Model parameters vector vector - [gamma_MC, gamma_X, E_Omega, E_X]
        params_calc_TE = [gamma_MC(width_index)/Consts.e_0, gamma_X_TE, abs(Fits{ii}.Dist_TE), Fits{ii}.TE];
        params_calc_TM = [gamma_MC(width_index)/Consts.e_0, gamma_X_TM, abs(Fits{ii}.Dist_TM), Fits{ii}.TM];
        
        % Calculating the anti-crossing diagrams
        [E_min_MCQW_TE_calc_abs,E_min_MCQW_TE_calc_r,E_min_MCQW_TE_calc_i,E_min_MCQW_TE_calc_amp] = CoupledOscillatorModelFunction(params_calc_TE, current_vec.Detuning.E_min_MC);
        [E_min_MCQW_TM_calc_abs,E_min_MCQW_TM_calc_r,E_min_MCQW_TM_calc_i,E_min_MCQW_TM_calc_amp] = CoupledOscillatorModelFunction(params_calc_TM, current_vec.Detuning.E_min_MC);
        
        % Saving results
        ACParams{ii}.CoupledOsc.Calc.width_index = width_index;
        ACParams{ii}.CoupledOsc.Calc.params_TE = params_calc_TE;
        ACParams{ii}.CoupledOsc.Calc.params_TM = params_calc_TM;
        N_TE = length(params_calc_TE(2:end))/3;
        ACParams{ii}.CoupledOsc.Calc.gamma_X_TE = params_calc_TE(2:N_TE+1); 
        ACParams{ii}.CoupledOsc.Calc.E_Omega_TE = params_calc_TE(N_TE+2:2*N_TE+1);
        ACParams{ii}.CoupledOsc.Calc.E_X_TE = params_calc_TE(2*N_TE+2:3*N_TE+1);
        N_TM = length(params_calc_TM(2:end))/3;
        ACParams{ii}.CoupledOsc.Calc.gamma_X_TM = params_calc_TM(2:N_TE+1); 
        ACParams{ii}.CoupledOsc.Calc.E_Omega_TM = params_calc_TM(N_TM+2:2*N_TM+1);
        ACParams{ii}.CoupledOsc.Calc.E_X_TM = params_calc_TM(2*N_TM+2:3*N_TM+1);
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TE_calc_abs = E_min_MCQW_TE_calc_abs;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TM_calc_abs = E_min_MCQW_TM_calc_abs;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TE_calc_r = E_min_MCQW_TE_calc_r;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TM_calc_r = E_min_MCQW_TM_calc_r;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TE_calc_i = E_min_MCQW_TE_calc_i;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TM_calc_i = E_min_MCQW_TM_calc_i;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TE_calc_amp = E_min_MCQW_TE_calc_amp;
        ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TM_calc_amp = E_min_MCQW_TM_calc_amp;
        ACParams{ii}.CoupledOsc.Calc.E_calc_vec = current_vec.Detuning.E_min_MC;
        ACParams{ii}.CoupledOsc.Calc.gamma_MC = gamma_MC;
    end
catch exc
    
end

%% Plotting

% try
%     disp(' - Plotting the results');
%     close all;
%     h_AC_TE = figure('Name', 'TE');
%     h_AC_TM = figure('Name', 'TM');
%     for (ii=1:eval(['length(MCParams' dim ' )']))
%         eval(['current_vec = MCParams' index ';']);
%         %current_vec.Detuning.E_min_MCQW_TE(current_vec.Detuning.E_min_MCQW_TE == 0) = nan;
%         %current_vec.Detuning.E_min_MCQW_TM(current_vec.Detuning.E_min_MCQW_TM == 0) = nan;
%         
%         figure(h_AC_TE);
%         subplot(1, length(MCParams), ii); hold on; box on;
%         plot(current_vec.Detuning.E_min_MC, current_vec.Detuning.E_min_MCQW_TE.', '.b', 'MarkerSize', 5); axis tight;
%         plot(current_vec.Detuning.E_min_MC, current_vec.Detuning.E_min_MC, ':r');
%         plot(ACParams{ii}.CoupledOsc.Fit.E_fit_vec, ACParams{ii}.CoupledOsc.Fit.E_min_MCQW_TE_fit_r, 'r');
%         %plot(ACParams{ii}.CoupledOsc.Calc.E_calc_vec, ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TE_calc_r, 'k');
%         if (~isempty(ACParams{ii}.TE))
%             plot(current_vec.Detuning.E_min_MC, repmat(ACParams{ii}.TE.', 1, length(current_vec.Detuning.E_min_MCQW_TE(1,:))), ':g', 'LineWidth', 1);
%             for (jj=1:length(ACParams{ii}.TE))
%                 text(current_vec.Detuning.E_min_MC(end), ACParams{ii}.TE(jj), [num2str(ACParams{ii}.TE(jj)) 'eV'], 'FontSize', 8);
%             end
%         end
%         if (m==1)
%             title([strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}'], 'FontSize', 8);
%         elseif (n==1)
%             title([strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}'], 'FontSize', 8);
%         end
%         if (ii~=1)
%             set(gca, 'YTickLabel', '');
%         end
%         if (ii==1)
%             ylabel('E [eV]');
%         end
%         xlabel('E [eV]');
%         axis([1.522 1.538 1.522 1.538]);
%         
%         figure(h_AC_TM);
%         subplot(1, length(MCParams), ii); hold on; box on;
%         plot(current_vec.Detuning.E_min_MC, current_vec.Detuning.E_min_MCQW_TM.', '.b', 'MarkerSize', 5); axis tight
%         plot(current_vec.Detuning.E_min_MC, current_vec.Detuning.E_min_MC, ':r');
%         plot(ACParams{ii}.CoupledOsc.Fit.E_fit_vec, ACParams{ii}.CoupledOsc.Fit.E_min_MCQW_TM_fit_r, 'r');
%         %plot(ACParams{ii}.CoupledOsc.Calc.E_calc_vec, ACParams{ii}.CoupledOsc.Calc.E_min_MCQW_TM_calc_r, 'k');
%         if (~isempty(ACParams{ii}.TM))
%             plot(current_vec.Detuning.E_min_MC, repmat(ACParams{ii}.TM.', 1, length(current_vec.Detuning.E_min_MCQW_TM(1,:))), ':g', 'LineWidth', 1);
%             for (jj=1:length(ACParams{ii}.TM))
%                 text(current_vec.Detuning.E_min_MC(end), ACParams{ii}.TM(jj), [num2str(ACParams{ii}.TM(jj)) 'eV'], 'FontSize', 8);
%             end
%         end
%         if (m==1)
%             title([strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}'], 'FontSize', 8);
%         elseif (n==1)
%             title([strrep(num2str(Params.gamma_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}'], 'FontSize', 8);
%         end
%         if (ii~=1)
%             set(gca, 'YTickLabel', '');
%         end
%         if (ii==1)
%             ylabel('E [eV]');
%         end
%         xlabel('E [eV]');
%         axis([1.522 1.538 1.522 1.538]);
%     end
% catch exc
%     
% end

% figure('Name', 'Differences');
% for (ii=1:eval(['length(MCParams' dim ' )']))
%     subplot(211); hold on; box on;
%     if (~isempty(ACParams{ii}.Diff_TE))
%         if (m==1)
%             semilogx(Params.N_DEG_vec(ii), abs(ACParams{ii}.Diff_TE)*1e3, '.b');
%         elseif (n==1)
%             semilogx(Params.gamma_vec(ii), abs(ACParams{ii}.Diff_TE)*1e3, '.b');
%         end
%         title('TE'); ylabel('\DeltaE [meV]');
%     end
%     subplot(212); hold on; box on;
%     if (~isempty(ACParams{ii}.Diff_TM))
%         if (m==1)
%             semilogx(Params.N_DEG_vec(ii), abs(ACParams{ii}.Diff_TM)*1e3, '.b');
%             xlabel('N_{DEG} [cm^{-2}]');
%         elseif (n==1)
%             semilogx(Params.gamma_vec(ii), abs(ACParams{ii}.Diff_TM)*1e3, '.b');
%             xlabel('\gamma [sec^{-1}]');
%         end
%         title('TM'); ylabel('\DeltaE [meV]');
%     end
% end

function F = ReshuffleVectors(E_min_MCQW, E_min_MC)

temp = E_min_MCQW;
temp(temp==0)=NaN;

[h,w]=size(temp);

flag = zeros(h,w);
new_temp = nan(h,w);

thr =6e-4;

for (jj=1:h)
    first_num_index = find(~isnan(temp(jj,:)),1,'first');
    %first_num_index = find(flag(jj,:)==0, 3, 'first');
    [i,j] = find(temp(:,first_num_index)~= new_temp(:,first_num_index));
    new_temp(jj,first_num_index)=temp(i(1),first_num_index);
    flag (jj,first_num_index)=1;
    kk=1;
    for (ii=first_num_index+1:w)
        
        %if (~isnan(new_temp(jj,ii-1)))
        dist_vec = abs(temp(:,ii)-new_temp(jj,ii-1));
        %else
        %    dist_vec = abs(temp(:,ii)-new_temp(jj,kk));
        %end
        
        if (isnan(min(dist_vec)))
            new_temp(jj,ii)=NaN;
            %flag (jj, ii)=1;
            continue;
        elseif (min(dist_vec)>=thr)
            new_temp(jj,ii)=NaN;
            kk=ii;
            %flag (jj, ii)=1;
            continue;
        end;
        
        [tt, index_min]=min(dist_vec);
        
        if (flag (index_min, ii)==0)
            new_temp(jj,ii)=temp(index_min, ii);
            flag (index_min, ii)=1;
        end;
    end;
end;

%new_temp(flag==0)=temp(flag==0)

%new_temp(new_temp==0)=NaN;

%new_temp = [new_temp;  resize(temp(flag==0), 1, length(new_temp(1,:)))];

flag(isnan(temp)) = nan;
new_temp(new_temp==0) = nan;

for (hh=1:h)
    if (length(flag(hh,flag(hh,:)==0))>20)
        vec = flag(hh, :);
        vec(vec==1)=nan;
        vec(vec==0)=temp(hh, flag(hh,:)==0);
        new_temp = [new_temp; vec];
    end
end

for (hh=1:length(new_temp(:,1)))
    nans(hh) = sum(isnan(new_temp(hh,:)));
end

new_temp = new_temp(nans<(3/4)*length(new_temp(1,:)), :);

for (hh=1:length(new_temp(:,1)))
    for (ww=2:w-1)
        if (~isnan(new_temp(hh,ww)) && isnan(new_temp(hh,ww-1)) && isnan(new_temp(hh,ww+1)))
            new_temp(hh,ww) = nan;
        end
    end
end

for (hh=1:length(new_temp(:,1)))
    first_non_nan = find(~isnan(new_temp(hh,:)), 1, 'first');
    last_non_nan = find(~isnan(new_temp(hh,:)), 1, 'last');
    new_temp(hh, 1:first_non_nan-1) = new_temp(hh,first_non_nan);
    new_temp(hh, last_non_nan+1:end) = new_temp(hh,last_non_nan);
end

% for (hh=1:length(new_temp(:,1)))
%     new_temp(hh,:) = interp1(E_min_MC, new_temp(hh,:), E_min_MC, 'pchip');
% end

[Y,I] = sort(new_temp(:,1), 'ascend');
F = new_temp(I,:);

function r = GetRandomNumber(min, max, x, y)

r = min + (max-min).*rand(x,y);

function gamma = GetMCLineWidth(MCParams)

global Consts;

%h = figure('Name','MC Reflection - Half Width');
r_MC_inv = 1-abs(MCParams.Detuning.r_MC);
E_vec = MCParams.Detuning.E_vec;
E_vec = min(E_vec):(max(E_vec)-min(E_vec))/20000:max(E_vec);

for (ii=1:length(r_MC_inv(:,1)))
    %current_vec = r_MC_inv(ii,:);
    current_vec = interp1(MCParams.Detuning.E_vec, r_MC_inv(ii,:), E_vec, 'pchip');
    [max_value, max_i] = max(current_vec);
    %    half_width_1 = interp1(current_vec, E_vec, max_value/2, 'pchip');
    %    if (half_width_1 > E_vec(max_i))
    %        half_width_2 = E_vec(max_i) - (half_width_1-E_vec(max_i));
    %    else
    %        half_width_2 = E_vec(max_i) + (E_vec(max_i)-half_width_1);
    %    end
    %    gamma(ii) = abs(half_width_1-half_width_2)/Consts.e_0;
    gamma(ii) = fwhm(E_vec, current_vec);
    
    %     figure(h); box on;
    %     plot(E_vec/Consts.e_0, current_vec, 'b'); hold on;
    %     plot(E_vec(max_i)/Consts.e_0, max_value, 'xr');
    %     plot((E_vec(max_i)+gamma(ii)/2)/Consts.e_0, max_value/2, 'xr', (E_vec(max_i)-gamma(ii)/2)/Consts.e_0, max_value/2, 'xr');
    %     %plot(half_width_1/Consts.e_0, max_value/2, 'xr', half_width_2/Consts.e_0, max_value/2, 'xr');
    %     hold off;
    %     xlabel('E [eV]'); ylabel('1-|r_{MC}|^2'); title(['\delta=' num2str(MCParams.Detuning.delta_vec(ii))]);
    %     hold on; axis tight;
end

function [gamma_vec_TE, gamma_vec_TM] = GetExcitonLineWidths(MC_Params, index)

global Consts;

h_TE = figure('Name', 'Line Width - TE');
h_TM = figure('Name', 'Line Width - TM');
num_peaks = str2double(input('-- How many peaks? ', 's'));

current_TE = 1-abs(MC_Params.Detuning.r_MCQW_TE(index,:));
current_TM = 1-abs(MC_Params.Detuning.r_MCQW_TM(index,:));
peaks_TE = interp1(MC_Params.Detuning.E_vec/Consts.e_0, current_TE, MC_Params.Detuning.E_min_MCQW_TE(:,index).', 'pchip');
peaks_TM = interp1(MC_Params.Detuning.E_vec/Consts.e_0, current_TM, MC_Params.Detuning.E_min_MCQW_TM(:,index).', 'pchip');

figure(h_TE); box on; hold on; axis([1.525, 1.54, 0, 0.3]);
plot(MC_Params.Detuning.E_vec/Consts.e_0, current_TE, 'b', MC_Params.Detuning.E_min_MCQW_TE(:,index).', peaks_TE, 'xr');
plot(MC_Params.Detuning.E_vec/Consts.e_0, repmat((peaks_TE/2).', 1, length(MC_Params.Detuning.E_vec)), ':g');
xlabel('E [eV]'); ylabel('1-|r_{TE}|');
for (ii=1:num_peaks)
    title(['TE - Mark ' num2str(ii) ' couple'])
    [x,y,but] = ginput(2);
    gamma_vec_TE(ii) = x(2)-x(1);
end

figure(h_TM); box on; hold on; axis([1.525, 1.54, 0, 0.3]);
plot(MC_Params.Detuning.E_vec/Consts.e_0, current_TM, 'b', MC_Params.Detuning.E_min_MCQW_TM(:,index).', peaks_TM, 'xr');
plot(MC_Params.Detuning.E_vec/Consts.e_0, repmat((peaks_TM/2).', 1, length(MC_Params.Detuning.E_vec)), ':g');
xlabel('E [eV]'); ylabel('1-|r_{TM}|');
for (ii=1:num_peaks)
    title(['TM - Mark ' num2str(ii) ' couple'])
    [x,y,but] = ginput(2);
    gamma_vec_TM(ii) = x(2)-x(1);
end

% but = 1; ACParams{ii}.TE = [];
% while (but == 1)
%     [x_TE,y_TE,but] = ginput(1);
%     if (but==1)
%         ACParams{ii}.TE = [ACParams{ii}.TE, y_TE];
%         plot(current_vec.Detuning.E_min_MC, y_TE*ones(1, length(current_vec.Detuning.E_min_MC)), ':r');
%         drawnow;
%     end
% end
% hold off;

