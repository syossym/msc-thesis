function Results = ReadDirectoryContents(dirname, Results, T)

global father_dirname;

dirname

filenames = dir(dirname);
jj = 1;

for (ii=1:length(filenames))
    if (filenames(ii).isdir && isempty(strfind(filenames(ii).name, '.')) && ~strcmp(filenames(ii).name, 'Backup'))
        if (~isempty(strfind(filenames(ii).name, 'K')))
            T = str2double(substr(filenames(ii).name, 0, length(filenames(ii).name)-1));
            eval(['Results.T' num2str(T) 'K.T = T;']);
        end
        
        Results = ReadDirectoryContents([dirname '\' filenames(ii).name], Results, T);
    else
        if (~isempty(strfind(filenames(ii).name, '.dat')))
            try
                filename_parts = strsplit('_', filenames(ii).name);
                file_type = filename_parts{1};
                if (strcmp(filename_parts{2}, 'lamp'))
                    spectral_offset = strrep(filename_parts{3}, '.dat', '');
                    file_type = 'Ref';
                    loaded_data = load([dirname '\' filenames(ii).name]);
                    %loaded_data(:,2) = loaded_data(:,2)./max(loaded_data(:,2));
                    eval(['Results.LampRef{jj}.Spectral_Offset = spectral_offset;']);
                    eval(['Results.LampRef{jj}.Data = loaded_data(60:910,:);']);                    
                else   
                    exposure = filename_parts{2};
                    spectral_offset = filename_parts{5};
                    delta_x = strrep(filename_parts{7}, '.dat', '');
                    
                    spectral_offset = strrep(spectral_offset, '-', 'm');
                    spectral_offset = strrep(spectral_offset, '+', 'p');
                    if (strcmp(spectral_offset, '0'))
                        spectral_offset = '000';
                    end
                    
                    if (strcmp(file_type, 'ref'))
                        file_type = 'Ref';
                    elseif (strcmp(file_type, 'pl'))
                        file_type = 'PL';
                    end
                    
                    eval(['Results.T' num2str(T) 'K.' file_type '.W' spectral_offset '.Exposure = exposure;']);
                    eval(['Results.T' num2str(T) 'K.' file_type '.W' spectral_offset '.Spectral_Offset = spectral_offset;']);
                    eval(['Results.T' num2str(T) 'K.' file_type '.W' spectral_offset '.dl(ii) = str2double(delta_x);']);
                    loaded_data = load([dirname '\' filenames(ii).name]);
                    %loaded_data(:,2) = loaded_data(:,2)./max(loaded_data(:,2));
                    eval(['Results.T' num2str(T) 'K.' file_type '.W' spectral_offset '.Data(ii,:) = loaded_data(60:910,2).'';']);
                    eval(['Results.T' num2str(T) 'K.' file_type '.W' spectral_offset '.E = loaded_data(60:910,1).'';']);
                end
                jj = jj + 1;
            catch exc
                disp(exc.message);
            end
        end
    end
end

