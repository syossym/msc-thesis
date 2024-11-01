function PlotExperimentalResults(Results)

global Consts;

%% Init

y_var.name = nan;
y_var.units = '';
y_var.num_type = 'reg';
colors = ['b', 'r', 'g', 'k', 'c', 'm', 'b', 'r', 'g', 'k', 'c', 'm'];

%% 300K

try
    if (IsStructField(Results.T300K, 'Ref') && IsStructField(Results.T300K, 'PL'))
        temp_num = 2;
    else
        temp_num = 1;
    end
    
    figure('Name', 'T=300K, Reflection, PL');
    if (IsStructField(Results.T300K, 'Ref'))
        fields = fieldnames(Results.T300K.Ref);
        index=length(fields):-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T300K.Ref.(fields{index(ii)}).E), max(Results.T300K.Ref.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T300K.Ref.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii);
            PlotWaterfall(Results.T300K.Ref.(fields{index(ii)}).E, Results.T300K.Ref.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            title(Results.T300K.Ref.(fields{index(ii)}).Spectral_Offset);
            %xlabel('E (eV)');
            if (ii==1)
                ylabel('Reflection (a.u.)');
            end
            if (temp_num ~= 1)
                set(gca, 'XTickLabel', '');
            else
                xlabel('E (4eV)');
            end
        end
    end
    if (IsStructField(Results.T300K, 'PL'))
        fields = fieldnames(Results.T300K.PL);
        index=length(fields):-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T300K.PL.(fields{index(ii)}).E), max(Results.T300K.PL.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T300K.PL.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii+length(fields));
            PlotWaterfall(Results.T300K.PL.(fields{index(ii)}).E, Results.T300K.PL.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            %title(Results.T300K.PL.(fields{index(ii)}).Spectral_Offset);
            xlabel('E (eV)');
            if (ii==1)
                ylabel('PL (a.u.)');
            end
        end
    end
    
    figure('Name', 'T=300K, Ref, PL Extrema');
    if (IsStructField(Results.T300K, 'Ref'))
        fields = fieldnames(Results.T300K.Ref);
        for(ii=1:length(fields))
            h = subplot(temp_num,length(fields),ii); box on; hold on;
            for (jj=1:length(Results.T300K.Ref.(fields{index(ii)}).Extrema))
                plot(Results.T300K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T300K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T300K.Ref.(fields{index(ii)}).Extrema{jj}, '.');
            end
            title(Results.T300K.Ref.(fields{index(ii)}).Spectral_Offset);
            %xlabel('\delta');
            if (ii==1)
                ylabel('E_{min}^{Ref} (eV)');
            end
            if (temp_num ~= 1)
                set(gca, 'XTickLabel', '');
            else
                xlabel('E (eV)');
            end
        end
    end
    if (IsStructField(Results.T300K, 'PL'))
        fields = fieldnames(Results.T300K.PL);
        for(ii=1:length(fields))
            h = subplot(temp_num,length(fields),ii+length(fields)); box on; hold on;
            for (jj=1:length(Results.T300K.PL.(fields{index(ii)}).Extrema))
                plot(Results.T300K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T300K.PL.(fields{index(ii)}).Extrema{jj})), Results.T300K.PL.(fields{index(ii)}).Extrema{jj}, '.');
            end
            title(Results.T300K.PL.(fields{index(ii)}).Spectral_Offset);
            xlabel('dl (a.u.)');
            if (ii==1)
                ylabel('E_{min}^{PL} (eV)');
            end
        end
    end
    
    %     figure('Name', 'T=300K, Ref, PL Extrema (joined)');
    %     for(ii=1:length(fields))
    %         h = subplot(1,length(fields),ii); box on; hold on;
    %         for (jj=1:length(Results.T300K.Ref.(fields{index(ii)}).Extrema))
    %             if (IsStructField(Results.T300K, 'Ref'))
    %                 plot(Results.T300K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T300K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T300K.Ref.(fields{index(ii)}).Extrema{jj}, '.b');
    %             end
    %             if (IsStructField(Results.T300K, 'PL'))
    %                 plot(Results.T300K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T300K.PL.(fields{index(ii)}).Extrema{jj})), Results.T300K.PL.(fields{index(ii)}).Extrema{jj}, '.r');
    %             end
    %         end
    %         axis tight;
    %         if (IsStructField(Results.T300K, 'Ref'))
    %             title(Results.T300K.Ref.(fields{index(ii)}).Spectral_Offset);
    %         end
    %         if (IsStructField(Results.T300K, 'PL'))
    %             title(Results.T300K.PL.(fields{index(ii)}).Spectral_Offset);
    %         end
    %         xlabel('Screw (a.u.)');
    %         if (ii==1)
    %             ylabel('E_{min} (eV)');
    %         end
    %     end
    
    figure('Name', 'T=300K, Ref, PL Extrema (joined^2)'); box on; hold on;
    title_text = [];
    fields = fieldnames(Results.T300K.Ref);
    for(ii=1:length(fields))
        for (jj=1:length(Results.T300K.Ref.(fields{index(ii)}).Extrema))
            if (IsStructField(Results.T300K, 'Ref'))
                plot(Results.T300K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T300K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T300K.Ref.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', '.', 'Color', colors(ii));
            end
        end
        axis tight;
        title_text = [Results.T300K.Ref.(fields{index(ii)}).Spectral_Offset, ', ' title_text];
        xlabel('Screw (a.u.)');
        if (ii==1)
            ylabel('E_{min} (eV)');
        end
    end
    fields = fieldnames(Results.T300K.PL);
    for(ii=1:length(fields))
        for (jj=1:length(Results.T300K.PL.(fields{index(ii)}).Extrema))
            if (IsStructField(Results.T300K, 'PL'))
                plot(Results.T300K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T300K.PL.(fields{index(ii)}).Extrema{jj})), Results.T300K.PL.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', 'x', 'Color', colors(ii));
            end
        end
        axis tight;
    end
    title(title_text);
    
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
    
    figure('Name', 'T=77K, Reflection, PL (1)');
    if (IsStructField(Results.T77K, 'Ref'))
        fields = fieldnames(Results.T77K.Ref);
        index=1:length(fields); %:-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T77K.Ref.(fields{index(ii)}).E), max(Results.T77K.Ref.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T77K.Ref.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii);
            PlotWaterfall(Results.T77K.Ref.(fields{index(ii)}).E, Results.T77K.Ref.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            title(Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset);
            xlabel('E (eV)');
            if (ii==1)
                ylabel('Reflection (a.u.)');
            end
            if (temp_num ~= 1)
                %set(gca, 'XTickLabel', '');
            else
                xlabel('E (4eV)');
            end
        end
    end
    if (IsStructField(Results.T77K, 'PL'))
        fields = fieldnames(Results.T77K.PL);
        index=1:length(fields); %:-1:1;
        %index=length(fields):-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T77K.PL.(fields{index(ii)}).E), max(Results.T77K.PL.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T77K.PL.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii+length(fields));
            PlotWaterfall(Results.T77K.PL.(fields{index(ii)}).E, Results.T77K.PL.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            title(Results.T77K.PL.(fields{index(ii)}).Spectral_Offset);
            xlabel('E (eV)');
            if (ii==1)
                ylabel('PL (a.u.)');
            end
        end
    end
    
%     figure('Name', 'T=77K, Reflection, PL (2)');
%     if (IsStructField(Results.T77K, 'Ref'))
%         fields = fieldnames(Results.T77K.Ref);
%         index=1:length(fields); %:-1:1;
%         for(ii=1:length(fields))
%             x_range = [min(Results.T77K.Ref.(fields{index(ii)}).E), max(Results.T77K.Ref.(fields{index(ii)}).E)];
%             [dl, dl_index] = sort(Results.T77K.Ref.(fields{index(ii)}).dl, 'descend');
%             y_var.value = dl;
%             
%             h = subplot(211); grid on; hold on;
%             %PlotWaterfall(Results.T77K.Ref.(fields{index(ii)}).E, Results.T77K.Ref.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
%             waterfall(Results.T77K.Ref.(fields{index(ii)}).E, 1:length(dl), Results.T77K.Ref.(fields{index(ii)}).Data(dl_index,:));
%             title(Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset);
%             xlabel('E (eV)');
%             if (ii==1)
%                 ylabel('Reflection (a.u.)');
%             end
%             if (temp_num ~= 1)
%                 %set(gca, 'XTickLabel', '');
%             else
%                 xlabel('E (4eV)');
%             end
%         end
%     end
%     if (IsStructField(Results.T77K, 'PL'))
%         fields = fieldnames(Results.T77K.PL);
%         index=1:length(fields); %:-1:1;
%         %index=length(fields):-1:1;
%         for(ii=1:length(fields))
%             x_range = [min(Results.T77K.PL.(fields{index(ii)}).E), max(Results.T77K.PL.(fields{index(ii)}).E)];
%             [dl, dl_index] = sort(Results.T77K.PL.(fields{index(ii)}).dl, 'descend');
%             y_var.value = dl;
%             
%             h = subplot(212); grid on; hold on;
%             %PlotWaterfall(Results.T77K.PL.(fields{index(ii)}).E, Results.T77K.PL.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
%             waterfall(Results.T77K.PL.(fields{index(ii)}).E, 1:length(dl), Results.T77K.PL.(fields{index(ii)}).Data(dl_index,:));
%             title(Results.T77K.PL.(fields{index(ii)}).Spectral_Offset);
%             xlabel('E (eV)');
%             if (ii==1)
%                 ylabel('PL (a.u.)');
%             end
%         end
%     end
    
    figure('Name', 'T=77K, Ref, PL Extrema');
    try
        if (IsStructField(Results.T77K, 'Ref'))
            fields = fieldnames(Results.T77K.Ref);
            for(ii=1:length(fields))
                h = subplot(temp_num,length(fields),ii); box on; hold on;
                for (jj=1:length(Results.T77K.Ref.(fields{index(ii)}).Extrema))
                    plot(Results.T77K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T77K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T77K.Ref.(fields{index(ii)}).Extrema{jj}, '.');
                end
                title(Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset);
                %xlabel('\delta');
                if (ii==1)
                    ylabel('E_{min}^{Ref} (eV)');
                end
                if (temp_num ~= 1)
                    set(gca, 'XTickLabel', '');
                else
                    xlabel('E (eV)');
                end
            end
        end
        if (IsStructField(Results.T77K, 'PL'))
            fields = fieldnames(Results.T77K.PL);
            for(ii=1:length(fields))
                h = subplot(temp_num,length(fields),ii+length(fields)); box on; hold on;
                for (jj=1:length(Results.T77K.PL.(fields{index(ii)}).Extrema))
                    plot(Results.T77K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T77K.PL.(fields{index(ii)}).Extrema{jj})), Results.T77K.PL.(fields{index(ii)}).Extrema{jj}, '.');
                end
                %title(Results.T77K.PL.(fields{index(ii)}).Spectral_Offset);
                xlabel('dl (a.u.)');
                if (ii==1)
                    ylabel('E_{min}^{PL} (eV)');
                end
            end
        end
    catch exc1
        disp(exc1.message);
    end
    
    %     figure('Name', 'T=77K, Ref, PL Extrema (joined)');
    %     for(ii=1:length(fields))
    %         h = subplot(1,length(fields),ii); box on; hold on;
    %         for (jj=1:length(Results.T77K.Ref.(fields{index(ii)}).Extrema))
    %             if (IsStructField(Results.T77K, 'Ref'))
    %                 plot(Results.T77K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T77K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T77K.Ref.(fields{index(ii)}).Extrema{jj}, '.b');
    %             end
    %             if (IsStructField(Results.T77K, 'PL'))
    %                 plot(Results.T77K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T77K.PL.(fields{index(ii)}).Extrema{jj})), Results.T77K.PL.(fields{index(ii)}).Extrema{jj}, '.r');
    %             end
    %         end
    %         axis tight;
    %         title(Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset);
    %         xlabel('Screw (a.u.)');
    %         if (ii==1)
    %             ylabel('E_{min} (eV)');
    %         end
    %     end
    
    figure('Name', 'T=77K, Ref, PL Extrema (joined^2)'); box on; hold on;
    title_text = [];
    fields = fieldnames(Results.T77K.Ref);
    try
        for(ii=1:length(fields))
            for (jj=1:length(Results.T77K.Ref.(fields{index(ii)}).Extrema))
                if (IsStructField(Results.T77K, 'Ref'))
                    plot(Results.T77K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T77K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T77K.Ref.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', '.', 'Color', colors(ii));
                end
            end
            axis tight;
            title_text = [Results.T77K.Ref.(fields{index(ii)}).Spectral_Offset, ', ' title_text];
            xlabel('Screw (a.u.)');
            if (ii==1)
                ylabel('E_{min} (eV)');
            end
        end
        fields = fieldnames(Results.T77K.PL);
        for(ii=1:length(fields))
            for (jj=1:length(Results.T77K.PL.(fields{index(ii)}).Extrema))
                if (IsStructField(Results.T77K, 'PL'))
                    plot(Results.T77K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T77K.PL.(fields{index(ii)}).Extrema{jj})), Results.T77K.PL.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', 'x', 'Color', colors(ii));
                end
            end
            axis tight;
        end
    catch exc1
        disp(exc1.message);
    end
    title(title_text);
    
catch exc2
    disp(exc2.message);
end

%% 2K

try
    if (IsStructField(Results.T2K, 'Ref') && IsStructField(Results.T2K, 'PL'))
        temp_num = 2;
    else
        temp_num = 1;
    end
    
    figure('Name', 'T=2K, Reflection, PL (1)');
    if (IsStructField(Results.T2K, 'Ref'))
        fields = fieldnames(Results.T2K.Ref);
        index=1:length(fields); %:-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T2K.Ref.(fields{index(ii)}).E), max(Results.T2K.Ref.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T2K.Ref.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii);
            PlotWaterfall(Results.T2K.Ref.(fields{index(ii)}).E, Results.T2K.Ref.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            title(Results.T2K.Ref.(fields{index(ii)}).Spectral_Offset);
            xlabel('E (eV)');
            if (ii==1)
                ylabel('Reflection (a.u.)');
            end
            if (temp_num ~= 1)
                %set(gca, 'XTickLabel', '');
            else
                xlabel('E (4eV)');
            end
        end
    end
    if (IsStructField(Results.T2K, 'PL'))
        fields = fieldnames(Results.T2K.PL);
        index=1:length(fields); %:-1:1;
        %index=length(fields):-1:1;
        for(ii=1:length(fields))
            x_range = [min(Results.T2K.PL.(fields{index(ii)}).E), max(Results.T2K.PL.(fields{index(ii)}).E)];
            [dl, dl_index] = sort(Results.T2K.PL.(fields{index(ii)}).dl, 'descend');
            y_var.value = dl;
            
            h = subplot(temp_num,length(fields),ii+length(fields));
            PlotWaterfall(Results.T2K.PL.(fields{index(ii)}).E, Results.T2K.PL.(fields{index(ii)}).Data(dl_index,:), h, 'b', x_range, y_var, 5);
            title(Results.T2K.PL.(fields{index(ii)}).Spectral_Offset);
            xlabel('E (eV)');
            if (ii==1)
                ylabel('PL (a.u.)');
            end
        end
    end
    
    figure('Name', 'T=2K, Ref, PL Extrema');
    try
        if (IsStructField(Results.T2K, 'Ref'))
            fields = fieldnames(Results.T2K.Ref);
            index=1:length(fields);
            for(ii=1:length(fields))
                h = subplot(temp_num,length(fields),ii); box on; hold on;
                for (jj=1:length(Results.T2K.Ref.(fields{index(ii)}).Extrema))
                    plot(Results.T2K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T2K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T2K.Ref.(fields{index(ii)}).Extrema{jj}, '.');
                end
                title(Results.T2K.Ref.(fields{index(ii)}).Spectral_Offset);
                %xlabel('\delta');
                if (ii==1)
                    ylabel('E_{min}^{Ref} (eV)');
                end
                if (temp_num ~= 1)
                    set(gca, 'XTickLabel', '');
                else
                    xlabel('E (eV)');
                end
            end
        end
        if (IsStructField(Results.T2K, 'PL'))
            fields = fieldnames(Results.T2K.PL);
            index=1:length(fields);
            for(ii=1:length(fields))
                h = subplot(temp_num,length(fields),ii+length(fields)); box on; hold on;
                for (jj=1:length(Results.T2K.PL.(fields{index(ii)}).Extrema))
                    plot(Results.T2K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T2K.PL.(fields{index(ii)}).Extrema{jj})), Results.T2K.PL.(fields{index(ii)}).Extrema{jj}, '.');
                end
                %title(Results.T2K.PL.(fields{index(ii)}).Spectral_Offset);
                xlabel('dl (a.u.)');
                if (ii==1)
                    ylabel('E_{min}^{PL} (eV)');
                end
            end
        end
    catch exc1
        disp(exc1.message);
    end
  
    figure('Name', 'T=2K, Ref, PL Extrema (joined^2)'); box on; hold on;
    title_text = [];
    fields = fieldnames(Results.T2K.Ref);
    index=1:length(fields);
    try
        for(ii=1:length(fields))
            for (jj=1:length(Results.T2K.Ref.(fields{index(ii)}).Extrema))
                if (IsStructField(Results.T2K, 'Ref'))
                    plot(Results.T2K.Ref.(fields{index(ii)}).dl(jj).*ones(1,length(Results.T2K.Ref.(fields{index(ii)}).Extrema{jj})), Results.T2K.Ref.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', '.', 'Color', colors(ii));
                end
            end
            axis tight;
            title_text = [Results.T2K.Ref.(fields{index(ii)}).Spectral_Offset, ', ' title_text];
            xlabel('Screw (a.u.)');
            if (ii==1)
                ylabel('E_{min} (eV)');
            end
        end
        fields = fieldnames(Results.T2K.PL);
        index=1:length(fields);
        for(ii=1:length(fields))
            for (jj=1:length(Results.T2K.PL.(fields{index(ii)}).Extrema))
                if (IsStructField(Results.T2K, 'PL'))
                    plot(Results.T2K.PL.(fields{index(ii)}).dl(jj).*ones(size(Results.T2K.PL.(fields{index(ii)}).Extrema{jj})), Results.T2K.PL.(fields{index(ii)}).Extrema{jj}, 'LineStyle', 'none', 'Marker', 'x', 'Color', colors(ii));
                end
            end
            axis tight;
        end
    catch exc1
        disp(exc1.message);
    end
    title(title_text);
    
catch exc2
    disp(exc2.message);
end

