clear all; close all; clc;

%% Constants (all MKS, except energy which is in eV)

global hbar e_0 a m_0;

hbar      = 1.054588757e-34; % Planck's constant over 2*pi
e_0       = 1.602189246e-19; % electronic charge
m_0       = 9.109534e-31;    % electron rest mass
c0        = 2.99792461e+8;   % speed of light
epsilon_0 = 8.854187827e-12; % permittivity of free space
h         = 6.626176583e-34; % Planck's constant
kb        = 1.380662e-23;    % Boltzmann constant
mu_b      = 9.274e-24;       % Bohr magneton


%% Luttinger-Kohn parameters

a     = 5.6605e-10;          % GaAs,AlAs lattice constant (meter)

%GaAs
g1  = 6.98; g2 = 2.06; g3 = 2.93;
w1  = (hbar^2/(2*m_0))*g1;
w2  = (hbar^2/(2*m_0))*g2;
w3  = (hbar^2/(2*m_0))*g3;
Ev  = 0;
Eg  = 1.426;                % band gap (eV)

% AlAs
g1  = 3.76; g2 = 0.82; g3 = 1.42;
a1  = (hbar^2/(2*m_0))*g1;
a2  = (hbar^2/(2*m_0))*g2;
a3  = (hbar^2/(2*m_0))*g3;

% Al(x)Ga(1-x)As
x   = 0.3;                % AlAs fraction in the alloy
b1  = ((1-x)*w1)+(x*a1);
b2  = ((1-x)*w2)+(x*a2);
b3  = ((1-x)*w3)+(x*a3);
% Evb = ((1-x)*0)+(x*0.75);
%Evb = ((1-x)*0)+(x*0.46);

%% Simulation

% z axis quantization
num_mono = 35;                   % number of monolayers in the well
W_mono   = 2.85;                 % monolayer width (A)
W        = 12; %num_mono*W_mono/10;   % well width (nm)
Wb       = 2*W;                  % cladding bulk width (nm)
Nw       = round(W*10);          % floor(1+W/(a*1e9));   % number of well layers
Nb       = 2*Nw;                 % number of bulk cladding layers
Np       = Nb+Nw+Nb;             % total number of layers

Z        = zeros(Np);            % z-axis

% Conduction band
% -----------------------------------------------------------
% Using the effective mass Sch. Eq. for single a particle (e)
% -----------------------------------------------------------
disp('-- Conduction band:');

hs_structure                        = [10*Wb, Nb, x ; 10*W, Nw, 0 ; 10*Wb, Nb, x];
[E_c,V_cb,V_vb,m_eff,con_profile]   = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
wf_c                                = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
z                                   = V_cb(:,1);
wf_c                                = NormaliseWfs(z, wf_c.');
m_eff_av                            = mean(m_eff(:,2));

dz   = z(2)-z(1);
Evb  = max(V_vb(:,2));

% Plotting
% ---------

figure(1); hold on;
%plot(z/1e-10, (Eg*e_0 + max(V_cb(:,2))*ones(length(z),1) - V_cb(:,2))/e_0, 'g');
plot(z/1e-10, (Eg*e_0 + V_cb(:,2))/e_0, 'g');
%V_vb = [zeros(1,Nb) Evb*ones(1,Nw) zeros(1,Nb)];
%plot(z/1e-10, V_vb - Evb*ones(size(V_vb)), 'g');
plot(z/1e-10, -V_vb(:,2)/e_0, 'g');

for(cb_index=1:length(E_c))
    plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, Eg + E_c(cb_index)*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), ':b');
    text(z(Nb+Nw+10)/1e-10, Eg + E_c(cb_index), sprintf('e%d',cb_index));
end
xlabel('z (Angstrom)'); ylabel('Energy (eV)');
title(['GaAs/Ga_{' num2str(x) '}Al_{' num2str(1-x) '}As, W_G_a_A_s=' num2str(W) 'nm, Eg(k=0)=' num2str(Eg) ...
    'eV: Energy levels in the structure']);
box on; grid on;

figure(2);
subplot(311);hold on;
for (cb_index=1:length(E_c))        % electrons
    plot(z/1e-10, abs(wf_c(cb_index,:)).^2);
    [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
    text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
end
ylabel('|F_e|^2');
title(['GaAs/Ga_{' num2str(x) '}Al_{' num2str(1-x) '}As, W_G_a_A_s=' num2str(W) ...
    'nm, Eg(k=0)=' num2str(Eg) 'eV: Probabity distributions of the conduction and valence bands']);
box on; grid on;

% Valence band
% ------------------------------
% Using Luttinger-Kohn approach
% ------------------------------
disp('-- Valence band:');

k_max = 0.02;   % (A^-1)
num_k = 5;     % size of k-axis

for (nk=1:(num_k+1))

    k(nk) = (nk-1)*(k_max/num_k);  % (A^-1)
    disp(['k=' num2str(k(nk)) '/' num2str(k_max) ' A^-1']);

    % Transverse k (x,y)
    l = 1; m = 0;
    lm  = sqrt((l^2)+(m^2));
    kx  = (l/lm)*k(nk)/1e-10;
    ky  = (m/lm)*k(nk)/1e-10;
    k2  = (kx^2)+(ky^2);
    k_t(nk) = sqrt(k2);

    % Upper part of the diagonalized Luttinger 2X2 Hamiltonian
    
    % H11
    H11_off   = (1/dz^2).*[(-b1+2*b2)*ones(1,Nb) , (-w1+2*w2)*ones(1,Nw-1) , (-b1+2*b2)*ones(1,Nb)];
    H11_diag2 = [0 H11_off] + [H11_off 0];
    temp_b    = (b1+b2)*k2+Evb;                  
    temp_w    = (w1+w2)*k2;
    temp_i    = 0.5*(temp_b+temp_w);
    H11_diag1 = [temp_b*ones(1,Nb) , temp_i , temp_w*ones(1,Nw-2) , temp_i , temp_b*ones(1,Nb)];
    
    H11       = diag(H11_diag1)+diag(H11_off,1)+diag(H11_off,-1)-diag(H11_diag2);
    
    % H22
    H22_off   = (1/dz^2).*[-(b1+2*b2)*ones(1,Nb) , -(w1+2*w2)*ones(1,Nw-1) , -(b1+2*b2)*ones(1,Nb)];
    H22_diag2 = [0 H22_off] + [H22_off 0];
    temp_b    = (b1-b2)*k2+Evb;                  
    temp_w    = (w1-w2)*k2;
    temp_i    = 0.5*(temp_b+temp_w);
    H22_diag1 = [temp_b*ones(1,Nb) , temp_i , temp_w*ones(1,Nw-2) , temp_i , temp_b*ones(1,Nb)];
    
    H22       = diag(H22_diag1)+diag(H22_off,1)+diag(H22_off,-1)-diag(H22_diag2);
    
    % H12
    H12_off   = -(1/(2*dz))*2*sqrt(3)*sqrt(k2).*[b3*ones(1,Nb) , w3*ones(1,Nw-1) , b3*ones(1,Nb)];
    temp_b    = 0.5*sqrt(3)*(b2+b3)*k2;
    temp_w    = 0.5*sqrt(3)*(w2+w3)*k2;
    temp_i    = 0.5*(temp_b+temp_w);
    H12_diag  = [temp_b*ones(1,Nb) , temp_i , temp_w*ones(1,Nw-2) , temp_i , temp_b*ones(1,Nb)];
    
    H12       = diag(H12_off,1)-diag(H12_off,-1)+diag(H12_diag);
    
    % H21
    H21_off   = (1/(2*dz))*2*sqrt(3)*sqrt(k2).*[b3*ones(1,Nb) , w3*ones(1,Nw-1) , b3*ones(1,Nb)];
    temp_b    = 0.5*sqrt(3)*(b2+b3)*k2;
    temp_w    = 0.5*sqrt(3)*(w2+w3)*k2;
    temp_i    = 0.5*(temp_b+temp_w);
    H21_diag  = [temp_b*ones(1,Nb) , temp_i , temp_w*ones(1,Nw-2) , temp_i , temp_b*ones(1,Nb)];
    
    H21       = diag(H21_off,1)-diag(H12_off,-1)+diag(H21_diag);
    
    H = -[H11, H12; H21, H22]./e_0;
        
    % The energies and states
    [V,D] = eig(H);
    D     = diag(D);
    [D,I] = sort(real(-D));
    D     = -(D)';
    V     = V(:,I);

    % Envelope functions
    F1{nk} = V(1:Np, 1:10);
    F2{nk} = V(Np+1:2*Np, 1:10);
    F_norm = sqrt(trapz(z,F1{nk}.^2+F2{nk}.^2,1));
    F1{nk} = F1{nk};
    F2{nk} = F2{nk};
    
    E_v_all(nk,:)  = D(1:10);
    %wf_v_all{nk}   = NormaliseWfs(z,V_f(1:Np,1:5));
    E_c_all(nk,:)  = ones(length(E_c),1).*Eg + E_c' + ...
        ones(length(E_c),1).*(((hbar^2)*k2)./(2*m_eff_av*e_0));

    % Plotting
    figure(2);
    
    subplot(312); hold on; box on; grid on;
    for (ii=1:length(F1{nk}(1,:)))       
        plot(z/1e-10, abs(F1{nk}(:,ii)).^2);
        [max_value,max_index] = max(abs(F1{nk}(:,ii)).^2);
        text(z(max_index)/1e-10, abs(F1{nk}(max_index,ii)).^2, sprintf('h%d',ii));
    end
    ylabel('|F^U_1|^2');
   
    subplot(313); hold on; box on; grid on;
    for (ii=1:length(F2{nk}(1,:)))       
        plot(z/1e-10, abs(F2{nk}(:,ii)).^2);
        [max_value,max_index] = max(abs(F2{nk}(:,ii)).^2);
        text(z(max_index)/1e-10, abs(F2{nk}(max_index,ii)).^2, sprintf('h%d',ii));
    end
    ylabel('|F^U_2|^2');
    xlabel('z (Angstrom)');
   
end

%clear V D;

%% Plotting

k     = k*10;     % (nm^-1)

% Potential profile
figure(1); hold on;
for (vb_index=1:length(E_v_all(1,:)))
    plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, E_v_all(1,vb_index).*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), ':r');
    text(z(Nb+Nw+10)/1e-10, E_v_all(1,vb_index), sprintf('h%d',vb_index));
end

% Dispersion
figure(3);
hold on;
for (ii=1:length(E_v_all(1,:)))
    plot(k_t, E_v_all(:,ii));
end
for (jj=1:length(E_c_all(1,:)))
    plot(k_t, E_c_all(:,jj));
end
title(['GaAs/Ga_{' num2str(x) '}Al_{' num2str(1-x) '}As, W=' num2str(W) 'nm: Conduction/Valence Bands E(k)']);
xlabel(['k_t : <' num2str(l) num2str(m) '>']);
ylabel('Energy (eV)')
grid on;


