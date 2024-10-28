function PlotFermiLevel(h_Band_Diagram, E_f, k_vec, T, E_c, E_h, N)

global Consts;

figure(h_Band_Diagram); hold on;
plot(k_vec, E_f*ones(length(k_vec))./Consts.e_0, 'g:');
text(min(k_vec), E_f./Consts.e_0, ['N=' num2str(N, '%10.0e') ' cm^{-2}: E_f=' num2str(E_f./Consts.e_0) ' eV'], 'FontSize', 8);

h_Fermi_Dirac = figure;
figure(h_Fermi_Dirac); 
subplot(211); hold on; box on;
for (ii=1:length(E_c))
    f_c = 1./(exp((E_c{ii}-E_f)/(Consts.k_B*T))+1);
    plot(k_vec, f_c);
    text(k_vec(1), f_c(1), ['cb=' num2str(ii)]); 
end
title('f_e');
xlabel('k_t [m^{-1}]');
subplot(212); hold on; box on;
for (ii=1:length(E_h))
    %f_h = 1./(exp((-E_h{ii}+E_f)/(Consts.k_B*T))+1);
    f_h = 1 - exp(-(E_h{ii}-E_f)/(Consts.k_B*T));
    plot(k_vec, f_h);
    text(k_vec(1), f_h(1), ['vb=' num2str(ii)]); 
end
title('f_h');
xlabel('k_t [m^{-1}]');

drawnow;

