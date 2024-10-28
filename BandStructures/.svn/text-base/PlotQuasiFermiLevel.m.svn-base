function PlotQuasiFermiLevel(h_Band_Diagram, F_c, F_h, k_vec, T, E_c, E_h, N, P)

global Consts;

k_c = F_c/(Consts.hbar*Consts.c);
k_h = F_h/(Consts.hbar*Consts.c);

figure(h_Band_Diagram); hold on;
subplot(222);
plot(k_vec, F_c*ones(length(k_vec))./Consts.e_0, 'g:');
text(min(k_vec), F_c./Consts.e_0, ['N=' num2str(N, '%10.0e') ' cm^{-2}: Ef_c=' num2str(F_c./Consts.e_0) ' eV'], 'FontSize', 8);
subplot(224);
plot(k_vec, F_h*ones(length(k_vec))./Consts.e_0, 'g:');
text(min(k_vec), F_h./Consts.e_0, ['P=' num2str(P, '%10.0e') ' cm^{-2}: Ef_h=' num2str(F_h./Consts.e_0) ' eV'], 'FontSize', 8);

h_Fermi_Dirac = figure;
figure(h_Fermi_Dirac); 
subplot(211); hold on; box on;
for (ii=1:length(E_c))
    f_c = 1./(exp((E_c{ii}-F_c)/(Consts.k_B*T))+1);
    plot(k_vec, f_c);
    text(k_vec(1), f_c(1), ['cb=' num2str(ii)]); 
end
%plot(k_c*ones(length(f_c)), linspace(0,0.01,length(f_c)), 'r:');
title(['N=' num2str(N, '%10.0e') ' cm^{-2}, T=' num2str(T) 'K, Ef_c=' num2str(F_c/Consts.e_0) ' eV']);
xlabel('k_t [m^{-1}]');
subplot(212); hold on; box on;
for (ii=1:length(E_h))
    f_h = 1./(exp((E_h{ii}-F_h)/(Consts.k_B*T))+1);
    plot(k_vec, f_h);
    text(k_vec(1), f_h(1), ['vb=' num2str(ii)]); 
end
%plot(k_h*ones(length(f_h)), linspace(0,0.01,length(f_h)), 'r:');
title(['P=' num2str(P, '%10.0e') ' cm^{-2}, T=' num2str(T) 'K, Ef_v=' num2str(F_h/Consts.e_0) ' eV']);
xlabel('k_t [m^{-1}]');

h_Fermi_Dirac_Factors = figure;
figure(h_Fermi_Dirac_Factors); 
subplot(211); hold on; box on;
for (ii=1:length(E_c))
    f_c = 1./(exp((E_c{ii}-F_c)/(Consts.k_B*T))+1);
    for (jj=1:length(E_h))
        f_h = 1./(exp((E_h{jj}-F_h)/(Consts.k_B*T))+1);
        factor_g = f_c+f_h-1; 
        plot(k_vec, factor_g);
        text(k_vec(1), factor_g(1), [num2str(ii) '-' num2str(jj)]); 
    end
end
title(['N=' num2str(N, '%10.0e') ' cm^{-2}, T=' num2str(T) 'K, E_f_c=' num2str(F_c/Consts.e_0) ' eV, E_f_v=' num2str(F_h/Consts.e_0) ' eV']);
subplot(212); hold on; box on;
for (ii=1:length(E_c))
    f_c = 1./(exp((E_c{ii}-F_c)/(Consts.k_B*T))+1);
    for (jj=1:length(E_h))
        f_h = 1./(exp((E_h{jj}-F_h)/(Consts.k_B*T))+1);
        factor_rsp = f_c.*f_h; 
        plot(k_vec, factor_rsp);
        text(k_vec(1), factor_rsp(1), [num2str(ii) '-' num2str(jj)]); 
    end
end
%plot(k_h*ones(length(f_h)), linspace(0,0.01,length(f_h)), 'r:');
xlabel('k_t [m^{-1}]');

drawnow;

