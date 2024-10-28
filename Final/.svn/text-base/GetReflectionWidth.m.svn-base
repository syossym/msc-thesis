function [Intensity,Widths] = GetReflectionWidth(MCParams, Params)

global Consts;

marker_types = ['.', 'o', 's', '^', 'v', 'd', '.', 'o', 's', '^', 'v', 'd'];
colors = ['b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm'];

close all; Fits = [];

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

options = optimset('TolFun',1e-5,'TolX',1e-5,'MaxIter',5e8,'MaxFunEvals',5e8, 'Display', 'off');
global E_0 Amp;

for (ii=1:eval(['length(MCParams' dim ' )']))
    eval(['current_mat = MCParams' index ';']);
    E_vec = current_mat.Detuning.E_vec./Consts.e_0;
    E_vec_interp = min(E_vec): (max(E_vec)-min(E_vec))/2000:max(E_vec);
    
    figure('Name', ['Linewidth & Intensity - N=' num2str(Params.N_DEG_vec(ii),'%1.0e')]);
    
    for (dd=1:length(current_mat.Detuning.delta_vec))
        
        current_t_TE = interp1(current_mat.Detuning.E_vec/Consts.e_0, smooth(1-abs(current_mat.Detuning.r_MCQW_TE(dd,:)),15), E_vec_interp, 'pchip');
        [maxs,mins] = peakdet(current_t_TE, 0.01, E_vec_interp);
        init_params = [maxs(:,1).', maxs(:,2).', 1e-4.*ones(1,length(maxs(:,1)))];
        %init_params = [1e-4.*length(maxs(:,1))];
        %E_0 = maxs(:,1).';
        %Amp = maxs(:,2).';
        
        %close all; 
        %figure(1); box on; hold on;
        %plot(E_vec_interp, current_t_TE, 'b');
        %axis([1.52,1.55,0,0.7]);
        
        fitted_params = lsqcurvefit(@LorentzianFitTargetFunction,init_params, E_vec_interp, current_t_TE, 0, Inf, options);
        
        %plot(E_vec_interp, LorentzianFitTargetFunction(fitted_params, E_vec_interp), 'r');
        %drawnow; 
        
        Widths.TE{ii,dd} = fitted_params(2*length(maxs(:,1))+1:end);
        Intensity.TE{ii,dd} = maxs(:,2).';
        
        current_t_TM = interp1(current_mat.Detuning.E_vec/Consts.e_0, smooth(1-abs(current_mat.Detuning.r_MCQW_TM(dd,:)),15), E_vec_interp, 'pchip');
        [maxs,mins] = peakdet(current_t_TM, 0.01, E_vec_interp);
        init_params = [maxs(:,1).', maxs(:,2).', 1e-4.*ones(1,length(maxs(:,1)))];
        
        Widths.TM{ii,dd} = fitted_params(2*length(maxs(:,1))+1:end);
        Intensity.TM{ii,dd} = maxs(:,2).';
        
        subplot(221); box on; hold on;
        for (jj=1:length(Intensity.TE{ii,dd}))
            plot(current_mat.Detuning.E_min_MC(dd), Intensity.TE{ii,dd}(jj),'Marker', '.', 'Color', colors(jj), 'LineStyle', 'none'); 
        end
        title('TE');
        ylabel('Transmission Peak');
        axis([min(current_mat.Detuning.E_min_MC), max(current_mat.Detuning.E_min_MC), 0, 1]);
        subplot(223); box on; hold on;
        for (jj=1:length(Widths.TE{ii,dd}))
            plot(current_mat.Detuning.E_min_MC(dd), 1e3*Widths.TE{ii,dd}(jj),'Marker', '.', 'Color', colors(jj), 'LineStyle', 'none'); 
        end
        ylabel('Linewidth (meV)'); xlabel('E_{CM} (eV)'); 
        axis([min(current_mat.Detuning.E_min_MC), max(current_mat.Detuning.E_min_MC), 0, 2]);
        
        subplot(222); box on; hold on;
        for (jj=1:length(Intensity.TM{ii,dd}))
            plot(current_mat.Detuning.E_min_MC(dd), Intensity.TM{ii,dd}(jj),'Marker', '.', 'Color', colors(jj), 'LineStyle', 'none'); 
        end
        title('TM');
        axis([min(current_mat.Detuning.E_min_MC), max(current_mat.Detuning.E_min_MC), 0, 1]);
        subplot(224); box on; hold on;
        for (jj=1:length(Widths.TM{ii,dd}))
            plot(current_mat.Detuning.E_min_MC(dd), 1e3*Widths.TM{ii,dd}(jj),'Marker', '.', 'Color', colors(jj), 'LineStyle', 'none'); 
        end
        xlabel('E_{CM} (eV)'); 
        axis([min(current_mat.Detuning.E_min_MC), max(current_mat.Detuning.E_min_MC), 0, 2]);
    end
%     ylabel('Linewidth (meV)');
%     xlabel('E_{MC} (eV)');
%     axis([min(current_mat.Detuning.E_min_MC), max(current_mat.Detuning.E_min_MC), 0, 2]);
    drawnow; 
end

