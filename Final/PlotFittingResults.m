function PlotFittingResults(Structure, Params, QWParams, MCParams, Fits)

marker_types = ['.', 'o', 's', '^', 'v', 'd'];
colors = ['b', 'r', 'g', 'k', 'c', 'm'];

figure(1);
for (ii=1:length(Params.N_DEG_vec))
    try
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.TE))
            plot(Params.N_DEG_vec(ii), Fits{ii}.TE(jj), '.');
        end
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_X_TE))
            plot(Params.N_DEG_vec(ii), Fits{ii}.CoupledOsc.Fit.E_X_TE(jj), 'ro');
        end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('E_{X}^{TE} (eV)');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.TM))
            plot(Params.N_DEG_vec(ii), Fits{ii}.TM(jj), '.');
        end
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_X_TM))
            plot(Params.N_DEG_vec(ii), Fits{ii}.CoupledOsc.Fit.E_X_TM(jj), 'ro');
        end
        title('(b)'); ylabel('E_{X}^{TM} (eV)'); xlabel('N_{2DEG} (cm^{-2})');
    catch exc
        continue;
    end
end

try
    fit_params_TE_mat = [];
    fit_params_TM_mat = [];
    for (ii=1:length(Params.N_DEG_vec))
        gamma_MC_TE(ii) = Fits{ii}.CoupledOsc.Fit.Params_TE(1);
        gamma_MC_TM(ii) = Fits{ii}.CoupledOsc.Fit.Params_TM(1);
        %fit_params_TE_mat = [fit_params_TE_mat ; Fits{ii}.CoupledOsc.Fit.Params_TE];
        %fit_params_TM_mat = [fit_params_TM_mat ; Fits{ii}.CoupledOsc.Fit.Params_TM];
    end
catch exc
    
end

figure(2);
semilogx(Params.N_DEG_vec(1:length(gamma_MC_TE)), gamma_MC_TE*1e3, 'b.', Params.N_DEG_vec(1:length(gamma_MC_TM)), gamma_MC_TM*1e3, 'r.'); 
hold on; plot(Params.N_DEG_vec, Fits{1}.CoupledOsc.Calc.params_TE(1)*1e3.*ones(1,length(Params.N_DEG_vec)), 'g');
xlabel('N_{2DEG} (cm^{-2})'); ylabel('$\gamma^{MC}$ (meV)', 'Interpreter', 'latex');

figure(3);
try
    for (ii=1:length(Params.N_DEG_vec))
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.Dist_TE))
            %plot(Params.N_DEG_vec(ii), Fits{ii}.Dist_TE(jj)*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        for (jj=1:length(Fits{ii}.Diff_TE))
            %plot(Params.N_DEG_vec(ii), abs(Fits{ii}.Diff_TE(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TE))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Fit.E_Omega_TE(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('$\hbar\Omega^{TE}$ (meV)', 'Interpreter', 'latex');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.Dist_TM))
            %plot(Params.N_DEG_vec(ii), Fits{ii}.Dist_TM(jj)*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        for (jj=1:length(Fits{ii}.Diff_TM))
            %plot(Params.N_DEG_vec(ii), abs(Fits{ii}.Diff_TM(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TM))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Fit.E_Omega_TM(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        title('(b)'); ylabel('$\hbar\Omega^{TM}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})');
    end
catch exc1
    % Do nothing
end

figure(4);
try
    for (ii=1:length(Params.N_DEG_vec))
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TE))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Calc.E_Omega_TE(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('$\hbar\Omega^{TE}$ (meV)', 'Interpreter', 'latex');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TM))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Calc.E_Omega_TM(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        title('(b)'); ylabel('$\hbar\Omega^{TM}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})');
    end
catch exc1
    % Do nothing
end

figure(5);
try
    for (ii=1:length(Params.N_DEG_vec))
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TE))
            plot(Params.N_DEG_vec(ii), (2*abs(Fits{ii}.CoupledOsc.Calc.E_Omega_TE(jj))*1e3).^2, 'r.');
        end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('$\hbar\Omega^{TE}$ (meV)', 'Interpreter', 'latex');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_Omega_TM))
            plot(Params.N_DEG_vec(ii), (2*abs(Fits{ii}.CoupledOsc.Calc.E_Omega_TM(jj))*1e3).^2, 'r.');
        end
        title('(b)'); ylabel('$(2\hbar\Omega^{TM})^2$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})');
    end
catch exc1
    % Do nothing
end

figure(6);
try
    for (ii=1:length(Params.N_DEG_vec))
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.TE))
            for (kk=1:length(Fits{ii}.TE))
                if (jj~=kk)
                    plot(Params.N_DEG_vec(ii), abs(Fits{ii}.TE(jj)-Fits{ii}.TE(kk))*1e3, '.');
                end
            end
        end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('$\Delta_{i,j}^{TE}$ (meV)', 'Interpreter', 'latex');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.TM))
            for (kk=1:length(Fits{ii}.TM))
                if (jj~=kk)
                    plot(Params.N_DEG_vec(ii), abs(Fits{ii}.TM(jj)-Fits{ii}.TM(kk))*1e3, '.');
                end
            end
        end
        title('(b)'); ylabel('$\Delta_{i,j}^{TM}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})');
    end
catch exc1
    % Do nothing
end

try
    for (ii=1:length(Params.N_DEG_vec))
        if (length(Fits{ii}.TE) >= 2)
            diff_E_X_12_Marked_TE(ii) =  abs(Fits{ii}.TE(1)-Fits{ii}.TE(2))*1e3;
        end
        if (length(Fits{ii}.TM) >= 2)
            diff_E_X_12_Marked_TM(ii) =  abs(Fits{ii}.TM(1)-Fits{ii}.TM(2))*1e3;
        end
    end
    
    for (ii=1:length(Params.N_DEG_vec))
        if (length(Fits{ii}.CoupledOsc.Fit.E_X_TE) >= 2)
            diff_E_X_12_Fitted_TE(ii) =  abs(Fits{ii}.CoupledOsc.Fit.E_X_TE(1)-Fits{ii}.CoupledOsc.Fit.E_X_TE(2))*1e3;
        end
        if (length(Fits{ii}.CoupledOsc.Fit.E_X_TM) >= 2)
            diff_E_X_12_Fitted_TM(ii) =  abs(Fits{ii}.CoupledOsc.Fit.E_X_TM(1)-Fits{ii}.CoupledOsc.Fit.E_X_TM(2))*1e3;
        end
    end
    
    figure(7);
    subplot(211); box on;
    plot(Params.N_DEG_vec(1:length(diff_E_X_12_Marked_TE)), diff_E_X_12_Marked_TE, '.-', Params.N_DEG_vec(1:length(diff_E_X_12_Marked_TM)), diff_E_X_12_Marked_TM, '.-r');
    ylabel('$\Delta_{e_1:hh_1,e_1:lh_1}^{Mark}$ (meV)', 'Interpreter', 'latex'); set(gca, 'XTickLabel', ''); title('(a)');
    subplot(212); box on;
    plot(Params.N_DEG_vec(1:length(diff_E_X_12_Fitted_TE)), diff_E_X_12_Fitted_TE, '.-', Params.N_DEG_vec(1:length(diff_E_X_12_Fitted_TM)), diff_E_X_12_Fitted_TM, '.-r');
    ylabel('$\Delta_{e_1:hh_1,e_1:lh_1}^{Fit}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})'); title('(b)');
catch exc1
    % Do nothing
end

try
    for (ii=1:length(Params.N_DEG_vec))
        if (length(Fits{ii}.CoupledOsc.Fit.E_X_TE) >= 2)
            diff_E_X_12_Fitted_TE(ii) =  abs(Fits{ii}.CoupledOsc.Fit.E_X_TE(1)-Fits{ii}.CoupledOsc.Fit.E_X_TE(2))*1e3;
        end
        if (length(Fits{ii}.CoupledOsc.Fit.E_X_TM) >= 2)
            diff_E_X_12_Fitted_TM(ii) =  abs(Fits{ii}.CoupledOsc.Fit.E_X_TM(1)-Fits{ii}.CoupledOsc.Fit.E_X_TM(2))*1e3;
        end
        
        if (length(Fits{ii}.CoupledOsc.Fit.E_Omega_TE) >= 1)
            Omega_E_X_12_Fitted_TE(ii) = abs(Fits{ii}.CoupledOsc.Fit.E_Omega_TE(1))*1e3;
        end
        if (length(Fits{ii}.CoupledOsc.Fit.E_Omega_TM) >= 1)
            Omega_E_X_12_Fitted_TM(ii) = abs(Fits{ii}.CoupledOsc.Fit.E_Omega_TM(1))*1e3;
        end
    end
    
    figure(8);
    subplot(211); box on;
    plot(Params.N_DEG_vec(1:length(diff_E_X_12_Fitted_TE)), diff_E_X_12_Fitted_TE, '.-', Params.N_DEG_vec(1:length(Omega_E_X_12_Fitted_TE)), 2.*Omega_E_X_12_Fitted_TE, '.-r');
    ylabel('$\Delta_{1,2}^{TE}, 2\hbar\Omega_{1,2}^{TE}$ (meV)', 'Interpreter', 'latex'); set(gca, 'XTickLabel', ''); title('(a)');
    subplot(212); box on;
    plot(Params.N_DEG_vec(1:length(diff_E_X_12_Fitted_TM)), diff_E_X_12_Fitted_TM, '.-', Params.N_DEG_vec(1:length(Omega_E_X_12_Fitted_TM)), 2.*Omega_E_X_12_Fitted_TM, '.-r');
    ylabel('$\Delta_{1,2}^{TM}, 2\hbar\Omega_{1,2}^{TM}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})'); title('(b)');
catch exc1
    % Do nothing
end

try
    for (ii=1:length(Params.N_DEG_vec))
        try
            Fits{ii}.CoupledOsc.Calc.params_TE;
        catch
            break;
        end
        gamma_MC_Fitted_TE(ii) = Fits{ii}.CoupledOsc.Fit.Params_TE(1);
        gamma_MC_Fitted_TM(ii) = Fits{ii}.CoupledOsc.Fit.Params_TM(1);
        %gamma_MC_Calc_TE(ii) = Fits{ii}.CoupledOsc.Calc.params_TE(1);
        %gamma_MC_Calc_TM(ii) = Fits{ii}.CoupledOsc.Calc.params_TM(1);
    end
    
    figure(9);
    subplot(211); box on;
    semilogx(Params.N_DEG_vec(1:length(gamma_MC_Fitted_TE)), gamma_MC_Fitted_TE*1e3, 'b.-', Params.N_DEG_vec(1:length(gamma_MC_Fitted_TM)), gamma_MC_Fitted_TM*1e3, 'r.-');
    ylabel('$\gamma_{MC}^{Fit}$ (meV)', 'Interpreter', 'latex'); set(gca, 'XTickLabel', ''); title('(a)');
    subplot(212); box on;
    semilogx(Params.N_DEG_vec(1:length(gamma_MC_Calc_TE)), gamma_MC_Calc_TE*1e3, 'b.-', Params.N_DEG_vec(1:length(gamma_MC_Calc_TM)), gamma_MC_Calc_TM*1e3, 'r.-');
    xlabel('N_{2DEG} (cm^{-2})'); ylabel('$\gamma_{MC}^{Calc}$ (meV)', 'Interpreter', 'latex'); title('(b)');
catch exc1
    % Do nothing
end

figure(10);
try
    for (ii=1:length(Params.N_DEG_vec))
        subplot(211); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.gamma_X_TE))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Fit.gamma_X_TE(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        %         for (jj=1:length(Fits{ii}.CoupledOsc.Calc.gamma_X_TE))
        %             plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Calc.gamma_X_TE(jj))*1e3, 'b.');
        %         end
        set(gca, 'XTickLabel', ''); title('(a)'); ylabel('$\gamma_{X}^{TE}$ (meV)', 'Interpreter', 'latex');
        subplot(212); box on; hold on;
        for (jj=1:length(Fits{ii}.CoupledOsc.Fit.gamma_X_TM))
            plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Fit.gamma_X_TM(jj))*1e3, 'Marker', marker_types(jj), 'Color', colors(jj));
        end
        %         for (jj=1:length(Fits{ii}.CoupledOsc.Calc.gamma_X_TM))
        %             plot(Params.N_DEG_vec(ii), abs(Fits{ii}.CoupledOsc.Calc.gamma_X_TM(jj))*1e3, 'b.');
        %         end
        title('(b)'); ylabel('$\gamma_{X}^{TM}$ (meV)', 'Interpreter', 'latex'); xlabel('N_{2DEG} (cm^{-2})');
    end
catch exc1
    % Do nothing
end

figure(11);
for (ii=1:length(Params.N_DEG_vec))
    subplot(5,round(length(Params.N_DEG_vec)/5),ii); box on; hold on;
    for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_min_MCQW_TE_fit_amp))
        plot(Fits{ii}.CoupledOsc.Fit.E_fit_vec, abs(Fits{ii}.CoupledOsc.Fit.E_min_MCQW_TE_fit_amp{jj}), 'b.', 'MarkerSize', 4.5);
    end
    %     for (jj=1:length(Fits{ii}.CoupledOsc.Fit.E_min_MCQW_TM_fit_amp))
    %         plot(Fits{ii}.CoupledOsc.Fit.E_fit_vec, abs(Fits{ii}.CoupledOsc.Fit.E_min_MCQW_TM_fit_amp{jj}), '.r', 'MarkerSize', 4.5);
    %     end
    if (ii==1)
        ylabel('$|\alpha_{P,i}|^2$', 'Interpreter', 'latex');
    end
    axis([min(Fits{ii}.CoupledOsc.Fit.E_fit_vec),max(Fits{ii}.CoupledOsc.Fit.E_fit_vec),0,1]);
    title([strrep(num2str(Params.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}'], 'FontSize', 8);
end