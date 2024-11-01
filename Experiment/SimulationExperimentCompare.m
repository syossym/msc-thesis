function SimulationExperimentCompare(Results, MCParams)

global Consts;

%% Init

y_var.name = nan;
y_var.units = '';
y_var.num_type = 'reg';
colors = ['b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm'];

%% 300K

try
    if (IsStructField(Results.T300K, 'Ref'))
        temp_num = 2;
    else
        temp_num = 1;
    end
    
    figure('Name', 'T=300K, Ref Extrema');
    h = subplot(211); box on; hold on; grid on;
    title_text = [];
    fields = fieldnames(Results.T300K.Ref);
    index=length(fields):-1:1;
    for(ii=1:length(fields))
        for (jj=1:length(Results.T300K.Ref.(fields{index(ii)}).Extrema))
            if (IsStructField(Results.T300K, 'Ref'))
                plot(2+Results.T300K.Ref.(fields{index(ii)}).dl(jj)*1.66.*ones(1,length(Results.T300K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T300K.Ref.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', '.', 'Color', colors(ii));
            end
        end
        axis tight;
        title_text = [Results.T300K.Ref.(fields{index(ii)}).Spectral_Offset, ', ' title_text];
        xlabel('Position (mm)');
        if (ii==1)
            ylabel('E_{min} (eV)');
        end
    end
    title(title_text);
    axis([0, 25, 1.5, 1.65]);
    
    subplot(212); box on; grid on; hold on;
    for (ii=length(MCParams.delta_vec):-1:1)
        plot(MCParams.delta_vec(ii), MCParams.Ref{ii}.mins(:,1), '.b');
        xlabel('\delta'); ylabel('E_{min} (eV)');
    end
    axis tight;
    axis([0.85, 1, 1.5, 1.65]);
    set(gca,'XDir','reverse')
    
catch exc1
    disp(exc1.message);
end

%% 77K

try
    if (IsStructField(Results.T77K, 'Ref') && IsStructField(Results.T77K, 'PL'))
        temp_num = 2;
    else
        temp_num = 1;
    end
    
    figure('Name', 'T=77K, Ref Extrema');
    subplot(211); box on; hold on; grid on;
    title_text = [];
    fields = fieldnames(Results.T77K.Ref);
    index=length(fields):-1:1;
    try
        for(ii=1:length(fields))
            for (jj=1:length(Results.T77K.Ref.(fields{index(ii)}).Extrema))
                if (IsStructField(Results.T77K, 'Ref'))
                    plot(2+Results.T77K.Ref.(fields{index(ii)}).dl(jj)*1.66.*ones(1,length(Results.T77K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T77K.Ref.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', '.', 'Color', colors(ii));
                end
            end
            axis tight;
            title_text = [Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset, ', ' title_text];
            xlabel('Position (mm)');
            if (ii==1)
                ylabel('E_{min} (eV)');
            end
        end
    catch exc1
        disp(exc1.message);
    end
    title(title_text);
    axis([0, 25, 1.5, 1.65]);
    
    subplot(212); box on; grid on; hold on;
    for (ii=length(MCParams.delta_vec):-1:1)
        plot(MCParams.delta_vec(ii), MCParams.Ref{ii}.mins(:,1), '.b');
        xlabel('\delta'); ylabel('E_{min} (eV)');
    end
    axis tight;
    axis([0.85, 1, 1.5, 1.65]);
    set(gca,'XDir','reverse')
catch exc2
    disp(exc2.message);
end
