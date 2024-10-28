function DataStruct = ExtractParameters(DataStruct, T, type)

figure('Name', type); box on;
for (ii=1:length(DataStruct.Data(:,1)))
    plot(DataStruct.E, DataStruct.Data(ii,:));
    xlabel('E [eV]'); title([type ' - Press the Return key to terminate the input']);
    hold on; axis tight;
    but = 1; DataStruct.Extrema{ii} = [];
    while (but == 1)
        [x,y,but] = ginput(1);
        if (but==1)
            DataStruct.Extrema{ii} = [DataStruct.Extrema{ii}, x];
            plot(x*ones(1,length(DataStruct.Data(ii,:))), DataStruct.Data(ii,:), ':r');
            drawnow;
        end
    end
    hold off;
end

close all;
