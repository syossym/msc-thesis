function Mapping = FitMCReflectionToExperiment(Params, delta_vec, Reflection, method)

global Consts;

options = optimset('TolFun',1e-15,'TolX',1e-15,'MaxIter',5e10,'MaxFunEvals',5e10, 'Display', 'on');
Mapping = [];

switch (method)
    
    case 'Calc',
        for (ii=1:length(Reflection))
            Params.T = Reflection{ii}.T;
            Structures{ii} = ReadSimulatedStructure(Params, Params.filename, Params.pathname);
        end
        
        for (ii=1:length(position))
            delta_current = GetDelta(Reflection, delta_vec, ii, Params, Structures);
            %delta_init = 0.9;
            %delta_fit = lsqcurvefit(@GetCavityMode, delta_init, position(ii), Reflection(ii), 0.9, 1, options);
            Mapping(ii,:) = [position(ii), delta_current]
        end
        
    case 'Interp',
        figure(1);
        for (jj=1:length(Reflection))
            position = Reflection{jj}.Cavity.Exp(:,1);
            for (ii=1:length(position))
                
                %                 if (jj==1)
                %                     E_exp_win = Reflection{jj}.Window.Exp(ii);
                %                 end
                E_exp_c = Reflection{jj}.Cavity.Exp(ii,2);
                
                % Window
                %                 if (jj==1)
                %                     if (E_exp_win <= max(Reflection{jj}.Window.Calc(:,1)))
                %                         delta_interp_win = interp1(Reflection{jj}.Window.Calc(:,1), Reflection{jj}.Window.Calc(:,2), E_exp_win, 'pchip');
                %                         Mapping.Window{jj}(ii,:) = [position(ii), delta_interp_win];
                %                         plot(Reflection{jj}.Window.Calc(:,1), Reflection{jj}.Window.Calc(:,2), '.', E_exp_win, delta_interp_win, 'xr');
                %                         drawnow;
                %                     else
                %
                %                         continue;
                %                     end
                %                 end
                
                % Cavity mode
                if (E_exp_c <= max(Reflection{jj}.Cavity.Calc(:,1)))
                    delta_interp_c = interp1(Reflection{jj}.Cavity.Calc(:,1), Reflection{jj}.Cavity.Calc(:,2), E_exp_c, 'pchip');
                    Mapping.Cavity{jj}(ii,:) = [position(ii), delta_interp_c];
                    plot(Reflection{jj}.Cavity.Calc(:,1), Reflection{jj}.Cavity.Calc(:,2), '.', E_exp_c, delta_interp_c, 'xr');
                    drawnow;
                else
                    
                    continue;
                end
            end
        end
end

function [cavity_mode_E, window_left_E] = GetCavityMode(delta, Params, Structure)

global Consts;

Params.delta = delta;
Ref = CalculateMCReflection(Structure, Params);
[maxs,mins] = peakdet(abs(Ref.r_MC), 0.0005, Params.E_exc/Consts.e_0);

[xmax,imax,xmin,imin] = extrema(diff(diff(mins(:,1))));
[m,i] = sort(xmax);
cavity_mode_E = mins(imax(1)+2,1);
window_left_E = mins(imax(1)+1,1);

figure(1);
plot(Ref.E_vec./Consts.e_0, abs(Ref.r_MC));
hold on;
plot(cavity_mode_E*ones(1,length(Ref.r_MC)), abs(Ref.r_MC), ':r', window_left_E*ones(1,length(Ref.r_MC)), abs(Ref.r_MC), ':g');
title(['\delta=' num2str(delta) ', T=' num2str(Params.T) 'K']);
hold off;
drawnow;
% if (length(mins(:,1))>0)
%     for (ii=1:length(mins(:,1)))
%         if (mins(ii,2) > 0.4 &&  mins(ii,2) < 0.5)
%             cavity_mode_E = mins(ii,1);
%             break;
%         else
%             cavity_mode_E = 0;
%         end
%     end
% else
%     cavity_mode_E = 0;
% end
%
% cavity_mode_E;

function delta = GetDelta(Reflection, delta_vec, pos_index, Params, Structures)

delta = zeros(1,length(Reflection));

for (jj=1:length(Reflection))
    Params.T = Reflection{jj}.T;
    for (ii=1:length(delta_vec))
        [cavity_mode_E, window_left_E] = GetCavityMode(delta_vec(ii), Params, Structures{jj});
        %         if (abs(Reflection{jj}.Cavity(pos_index)-cavity_mode_E)^2 < 1e-6)
        %             delta(jj) = delta_vec(ii);
        %             disp('break');
        %             break;
        %         end
        Reflection{jj}.Window(pos_index)
        window_left_E
        if (abs(Reflection{jj}.Window(pos_index)-window_left_E)^2 < 1e-6)
            delta(jj) = delta_vec(ii);
            disp('break');
            break;
        end
    end
end