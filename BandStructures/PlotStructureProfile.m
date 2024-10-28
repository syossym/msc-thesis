function h_structure_profile = PlotStructureProfile(parameters_struct, T, save)

global Consts;

%% Figures init

h_structure_profile = figure('Name','Structure Profile');

%% Plotting

E_g = GetMaterialBandGap(parameters_struct.Materials{1},T);
for (ii=1:length(length(parameters_struct.Materials)))
    if (parameters_struct.Materials{ii}.E_g<E_g)
       E_g = GetMaterialBandGap(parameters_struct.Materials{ii},T);
    end
end
separator_vec = -max(parameters_struct.v_h_profile)/Consts.e_0: 0.01 : max(parameters_struct.v_e_profile)/Consts.e_0+E_g;
temp = 0;
figure(h_structure_profile); box on; hold on;
plot(parameters_struct.z_grid/1e-10, zeros(1,length(parameters_struct.z_grid)), 'r:');
plot(parameters_struct.z_grid/1e-10, E_g.*ones(1,length(parameters_struct.z_grid)), 'r:');
text(0, E_g, ['E_g=' num2str(E_g) 'eV'], 'FontSize', 8);
plot(parameters_struct.z_grid/1e-10, E_g + parameters_struct.v_e_profile./Consts.e_0, 'b', ...
     parameters_struct.z_grid/1e-10, -parameters_struct.v_h_profile./Consts.e_0, 'b');
 
for (ii=1:length(parameters_struct.Materials))
   text(temp + parameters_struct.Materials{ii}.Width/2,E_g/2,[parameters_struct.Materials{ii}.Name ': x=' num2str(parameters_struct.Materials{ii}.x) ], 'FontSize', 10);
   temp = temp + parameters_struct.Materials{ii}.Width;
   plot(temp*ones(1,length(separator_vec)), separator_vec, 'g:');
end
xlabel('z [Angstrom]');
ylabel('E [eV]');
title(['Simulated Structure - T=' num2str(T) 'K']);

%% Optional figure saving

if (save)
    dirPath = uigetdir('c:\', 'Choose Directory');
    if (dirPath ~= 0)
        saveas(h_structure_profile, [dirPath '\StructureProfile'], 'fig');
        exportfig(h_structure_profile, [dirPath '\StructureProfile'], 'FontMode', 'scaled', 'FontSize', 1, 'color', 'cmyk', 'Format', 'eps');
    end
end