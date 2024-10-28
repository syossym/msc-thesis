function PlotOpticalSimulationParams(QStruct, OpticalParams, FullStructureParams, params)
%
% This function plots the optical simulation results.
%
%   Input: 'QStruct' - the simulated quantum structure.
%          'OpticalParams' - the optical susceptability (and derivatives)
%                            calculation results.
%          'FullStructureParams' - the reflection spectra simulation
%                                  results.
%          'params' - simulaion parameters.
%
%   Output: -
%
% Tested: Matlab 7.8.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

if (length(QStruct.N_DEG_vec) > 1 && length(params.gamma) == 1)
    sim_type = 1;
elseif (length(QStruct.N_DEG_vec) == 1 && length(params.gamma) > 1)
    sim_type = 2;
end

colors_vec = ['b', 'r', 'g', 'k'];

switch (sim_type)

    case 1,    % multiple DEG concentrations and single gamma value

        h_AntiCrossing_all = figure('Name', 'Anticrossing Diagrams');
        h_AntiCrossing_separate = figure('Name', 'Anticrossing Diagrams');
        h_Susceptability = figure('Name', 'Optical Susceptability');

        for (ii=1:length(FullStructureParams))
            arranged_curves = ArrangeAntiCrossingData(FullStructureParams{ii}.E_min_r_MC,FullStructureParams{ii}.E_min_r_MCQW);
            for (jj=1:length(arranged_curves))
                if (~isempty(arranged_curves{jj}))
                    figure(h_AntiCrossing_all); hold on; box on;
                    plot(arranged_curves{jj}.E_MC, arranged_curves{jj}.E_MCQW, 'b.', 'MarkerSize', 4); %);
                    plot(arranged_curves{jj}.E_MC, arranged_curves{jj}.E_MC, 'k');
                    xlabel('E_{MC} [eV]'); ylabel('E [eV]');
                    title(['\gamma=' strrep(num2str(params.gamma,'%1.0e'), 'e+0', 'x10^{') '}sec^{-1}']);
                    axis([1.52,1.54,1.52,1.54]);

                    figure(h_AntiCrossing_separate);
                    subplot(ceil(length(FullStructureParams)/3), 3, ii); hold on; box on;
                    plot(arranged_curves{jj}.E_MC, arranged_curves{jj}.E_MCQW, 'b.', 'MarkerSize', 4); %[colors_vec(jj) '.']);
                    plot(arranged_curves{jj}.E_MC, arranged_curves{jj}.E_MC, 'r');
                    %xlabel('E_{MC} [eV]'); ylabel('E [eV]');
                    title(['N_{2DEG}=' strrep(num2str(QStruct.N_DEG_vec(ii),'%1.0e'), 'e+0', 'x10^{') '}cm^{-2}']);
                    axis([1.52,1.54,1.52,1.54]);
                end
            end

            figure(h_Susceptability);
            
        end

    case 2,    % single DEG concentration and multiple gamma values

end

function f = SortMatrix(f)

[len_x,len_y] = size(f);
if (len_x > len_y)
    f = f.';
end

for (ii=1:len_y)
    f(:,ii) = sort(f(:,ii));
end

function f = ArrangeAntiCrossingData(E_r_m_MC,E_r_m_MCQW)

% Remove the outliers
for (ii=1:length(E_r_m_MCQW(:,1)))
    f{ii}.E_MC = E_r_m_MC(E_r_m_MCQW(ii,:)>1.5);
    f{ii}.E_MCQW = E_r_m_MCQW(ii, E_r_m_MCQW(ii,:)>1.5);
    f{ii}.Max = max(f{ii}.E_MCQW);
    f{ii}.Min = min(f{ii}.E_MCQW);
end




