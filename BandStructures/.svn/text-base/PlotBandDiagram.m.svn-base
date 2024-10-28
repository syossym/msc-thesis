function h_band_diagram = PlotBandDiagram(k_vec, E_c, E_v, x, T, L_z, well_material, save)

global Consts;

%% Figures init

h_band_diagram = figure('Name','Band Diagram');

%% Figures pre-processing

if (E_v{1}(1) > 0)
   for (ii=1:length(E_v))
      E_v{ii} = -E_v{ii}; 
   end
end
E_g = GetMaterialBandGap(well_material,T);

%% Plotting

max_index = ceil(length(k_vec)/2);

figure(h_band_diagram); 
subplot(2,2, [1 3]);
hold on;
axis([0, k_vec(max_index), -max(abs(E_v{end}(1:max_index)))/Consts.e_0, max(E_c{end}(1:max_index))/Consts.e_0]);
for (ii=1:length(E_c))
    plot(k_vec, E_c{ii}/Consts.e_0);
end
for (ii=1:length(E_v))
    plot(k_vec, E_v{ii}/Consts.e_0);
end
plot(k_vec, E_g.*ones(length(k_vec)), 'r:', k_vec, zeros(length(k_vec)), 'r:');
text(min(k_vec), E_g, ['E_g=' num2str(E_g) 'eV'], 'FontSize', 8);
title(['GaAs/Al_{' num2str(x) '}Ga_{' num2str(1-x) '}As, L_z=' num2str(L_z/1e-10) 'A, T=' num2str(T) 'K']);
xlabel('k_t [m^{-1}]');
ylabel('E [eV]');
box on;
subplot(222);
hold on;
axis([0, k_vec(max_index), E_g, max(E_c{end}(1:max_index))/Consts.e_0]);
for (ii=1:length(E_c))
    plot(k_vec, E_c{ii}/Consts.e_0);
end
text(min(k_vec), E_g, ['E_g=' num2str(E_g) 'eV'], 'FontSize', 8);
xlabel('k_t [m^{-1}]');
ylabel('E [eV]');
title('Conduction subbands');
box on;
subplot(224);
hold on;
axis([0, k_vec(max_index), -max(abs(E_v{end}(1:max_index))/Consts.e_0), 0]);
for (ii=1:length(E_v))
    plot(k_vec, E_v{ii}/Consts.e_0);
end
xlabel('k_t [m^{-1}]');
ylabel('E [eV]');
title('Valence subbands');
box on;

drawnow;

%% Optional figure saving

if (save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_band_diagram, [dirPath '\BandDiagram'], 'fig');
        exportfig(h_band_diagram, [dirPath '\BandDiagram'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
    end
end