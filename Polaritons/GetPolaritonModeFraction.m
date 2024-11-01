function ModeFraction = GetPolaritonModeFraction(Params, MCParams, Fits)

global Consts;

max_polariton_num = 7;

for (ii=1:length(Fits))
    current_Fit = Fits{ii};
    current_MC_Reflection = MCParams{ii,1};
    
    current_vec_TE = ReshuffleVectors(current_MC_Reflection.Detuning.E_min_MCQW_TE, current_MC_Reflection.Detuning.E_min_MC);
    current_vec_TM = ReshuffleVectors(current_MC_Reflection.Detuning.E_min_MCQW_TM, current_MC_Reflection.Detuning.E_min_MC);
    current_vec_TE(isnan(current_vec_TE)) = 0;
    current_vec_TM(isnan(current_vec_TM)) = 0;
    E_min_MC = sort(current_MC_Reflection.Detuning.E_min_MC, 'ascend');
    
    try
        current_Fit.CoupledOsc.Fit;
    catch exc
        continue;
    end
    
    [D_abs_TE, D_real_TE, D_imag_TE, V_vec_TE] = CoupledOscillatorModelFunction(current_Fit.CoupledOsc.Fit.Params_TE, current_Fit.CoupledOsc.Fit.E_fit_vec);
    [D_abs_TM, D_real_TM, D_imag_TM, V_vec_TM] = CoupledOscillatorModelFunction(current_Fit.CoupledOsc.Fit.Params_TM, current_Fit.CoupledOsc.Fit.E_fit_vec);
    
    for (gg=1:length(current_MC_Reflection.Detuning.r_MCQW_TE(:,1)))
        polariton_linewidth_TE{gg} = multfwhm(current_MC_Reflection.Detuning.E_vec/Consts.e_0,(1 - abs(current_MC_Reflection.Detuning.r_MCQW_TE(gg,:))));
        polariton_linewidth_TM{gg} = multfwhm(current_MC_Reflection.Detuning.E_vec/Consts.e_0,(1 - abs(current_MC_Reflection.Detuning.r_MCQW_TM(gg,:))));
    end
    
    %     figure(1);
    %     for (pp=1:length(V_vec_TE))
    %         subplot(length(Fits),4,(ii-1)*4+pp); box on;
    %         plot(current_Fit.CoupledOsc.Fit.E_fit_vec, abs(V_vec_TE{pp}.').^2);
    %         set(gca, 'XTickLabel', '');
    %         if (pp~=1), set(gca, 'YTickLabel', ''); end
    %         axis([1.522 1.538 0 1]);
    %     end
    %     figure(2);
    %     for (pp=1:length(V_vec_TM))
    %         subplot(length(Fits),4,(ii-1)*4+pp); box on;
    %         plot(current_Fit.CoupledOsc.Fit.E_fit_vec, abs(V_vec_TM{pp}.').^2);
    %         set(gca, 'XTickLabel', '');
    %         if (pp~=1), set(gca, 'YTickLabel', ''); end
    %         axis([1.522 1.538 0 1]);
    %     end
    
    % Plotting the anti-crossing and the polariton fractions
    figure('Name', ['Polariton Fraction - ' num2str(Params.N_DEG_vec(ii),'%1.0e')]);
    TE_indices = [1,3,5,7];
    TM_indices = [9,11,13,15];
    subplot(8,2,TE_indices); box on;
    plot(current_MC_Reflection.Detuning.E_min_MC, current_vec_TE.', '.',...
        current_MC_Reflection.Detuning.E_min_MC, current_MC_Reflection.Detuning.E_min_MC, ':r'); axis tight;
    title('TE'); ylabel('E (eV)');
    set(gca, 'XTickLabel', '');
    for (pp=1:length(V_vec_TE))
        if (pp<=4)
            subplot(8,2, TE_indices(pp)+1); box on;
            plot(current_Fit.CoupledOsc.Fit.E_fit_vec, abs(V_vec_TE{pp}.').^2);
            set(gca, 'XTickLabel', '');
            axis([1.522 1.538 0 1]);
            ylabel(['|\alpha_{' num2str(pp) ',i}|^2']);
        end
    end
    subplot(8,2,TM_indices); box on;
    plot(current_MC_Reflection.Detuning.E_min_MC, current_vec_TM.', '.',...
        current_MC_Reflection.Detuning.E_min_MC, current_MC_Reflection.Detuning.E_min_MC, ':r'); axis tight;
    title('TM'); xlabel('E_{CM} (eV)'); ylabel('E (eV)');
    for (pp=1:length(V_vec_TM))
        if (pp<=4)
            subplot(8,2, TM_indices(pp)+1); box on;
            plot(current_Fit.CoupledOsc.Fit.E_fit_vec, abs(V_vec_TM{pp}.').^2);
            if (pp~=length(V_vec_TM))
                set(gca, 'XTickLabel', '');
            else
                xlabel('E (eV)');
            end
            axis([1.522 1.538 0 1]);
            ylabel(['|\alpha_{' num2str(pp) ',i}|^2']);
        end
    end
    
    figure('Name', ['Polariton Fraction (TM) - ' num2str(Params.N_DEG_vec(ii),'%1.0e')]);
    TM_indices = [1,3,5,7,9,11,13,15];
    subplot(8,2,TM_indices); box on;
    plot(current_MC_Reflection.Detuning.E_min_MC, current_vec_TM.', '.',...
        current_MC_Reflection.Detuning.E_min_MC, current_MC_Reflection.Detuning.E_min_MC, ':r'); axis tight;
    title('TM'); xlabel('E_{CM} (eV)'); ylabel('E (eV)');
    for (pp=1:length(V_vec_TM))
        if (pp<=8)
            subplot(8,2, TM_indices(pp)+1); box on;
            plot(current_Fit.CoupledOsc.Fit.E_fit_vec, abs(V_vec_TM{pp}.').^2);
            if (pp~=length(V_vec_TM))
                set(gca, 'XTickLabel', '');
            else
                xlabel('E (eV)');
            end
            axis([1.52 1.545 0 1]);
            ylabel(['|\alpha_{' num2str(pp) ',i}|^2']);
        end
    end
    
    figure('Name', ['Imag - ' num2str(Params.N_DEG_vec(ii),'%1.0e')]);
    subplot(2,2,1); box on;
    plot(current_MC_Reflection.Detuning.E_min_MC, current_vec_TE.', '.',...
        current_MC_Reflection.Detuning.E_min_MC, current_MC_Reflection.Detuning.E_min_MC, ':r'); axis tight;
    title('TE'); ylabel('E (eV)');
    set(gca, 'XTickLabel', '');
    subplot(222); box on;
    plot(current_Fit.CoupledOsc.Fit.E_fit_vec(end:-1:1), -1e3*(D_imag_TE(:, end:-1:1)).'); axis tight;
    %     for (gg=1:length(polariton_linewidth_TE))
    %         hold on; plot(E_min_MC(gg), polariton_linewidth_TE{gg}*1e3, '.b');
    %     end
    set(gca, 'XTickLabel', ''); ylabel('Linewidth (TE) (meV)');
    subplot(2,2,3); box on;
    plot(current_MC_Reflection.Detuning.E_min_MC, current_vec_TE.', '.',...
        current_MC_Reflection.Detuning.E_min_MC, current_MC_Reflection.Detuning.E_min_MC, ':r'); axis tight;
    title('TM'); xlabel('E_{CM} (eV)'); ylabel('E (eV)');
    subplot(224); box on;
    plot(current_Fit.CoupledOsc.Fit.E_fit_vec(end:-1:1), -1e3*(D_imag_TM(:, end:-1:1)).'); axis tight;
    %     for (gg=1:length(polariton_linewidth_TM))
    %         hold on; plot(E_min_CM(gg), polariton_linewidth_TM{gg}*1e3, '.b');
    %     end
    ylabel('Linewidth (TM) (meV)'); xlabel('E_{CM} (meV)');
end

ModeFraction = [];