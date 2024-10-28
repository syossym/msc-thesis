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
C     = (hbar^2)/(2*m_0*e_0*(a^2));

%GaAs
g1_w  = 6.98; g2_w= 2.06; g3_w = 2.93;
w1  = (hbar^2)*g1_w/(2*m_0*e_0*(a^2));
w2  = (hbar^2)*g2_w/(2*m_0*e_0*(a^2));
w3  = (hbar^2)*g3_w/(2*m_0*e_0*(a^2));
Ev  = 0;
Eg  = 1.426;                % band gap (eV)

% AlAs
g1_b  = 3.76; g2_b = 0.82; g3_b = 1.42;
a1  = (hbar^2)*g1_b/(2*m_0*e_0*(a^2));
a2  = (hbar^2)*g2_b/(2*m_0*e_0*(a^2));
a3  = (hbar^2)*g3_b/(2*m_0*e_0*(a^2));

% Al(x)Ga(1-x)As
x   = 0.25;                % AlAs fraction in the alloy
b1  = ((1-x)*w1)+(x*a1);
b2  = ((1-x)*w2)+(x*a2);
b3  = ((1-x)*w3)+(x*a3);
% Evb = ((1-x)*0)+(x*0.75);
%Evb = ((1-x)*0)+(x*0.46);

%% Simulation

% z-direction quantization
num_mono = 23;                   % number of monolayers in the well
W_mono   = 2.85;                 % monolayer width (A)
W        = num_mono*W_mono/10;   % well width (nm)
Wb       = 2*W;                  % cladding bulk width (nm)
Nw       = round(W*10);          % floor(1+W/(a*1e9));   % number of well layers
Nb       = 2*Nw;                 % number of bulk cladding layers
Np       = Nb+Nw+Nb;             % total number of layers

Z        = zeros(Np);            % z-axis

% Luttinger parameter vectors
g1_vec = [g1_b*ones(1,Nb), g1_w*ones(1,Nw), g1_b*ones(1,Nb)];
g2_vec = [g2_b*ones(1,Nb), g2_w*ones(1,Nw), g2_b*ones(1,Nb)];
g3_vec = [g3_b*ones(1,Nb), g3_w*ones(1,Nw), g3_b*ones(1,Nb)];

% Conduction band - using the effective mass Sch. eq. for single particle (e)
disp('-- Conduction band:');

%hs_structure = [10*(Nb-1)*(a*1e9), Nb, x ; 10*W, Nw, 0 ; 10*(Nb-1)*(a*1e9), Nb, x];   % the width are given in A
hs_structure = [10*Wb, Nb, x ; 10*W, Nw, 0 ; 10*Wb, Nb, x];
[E_c,V_cb,V_vb,m_eff,con_profile] = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
wf_c = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
m_eff_av = mean(m_eff(:,2));
z = V_cb(:,1);
dz  = z(2)-z(1);
Evb = max(V_vb(:,2))/e_0;

% Plotting

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

% Valence band - using Luttinger-Kohn approach
disp('-- Valence band:');

k_max = 0.1;   % (A^-1)
num_k = 50;     % size of k-axis

for (nk=1:(num_k+1))

    k(nk) = (nk-1)*(k_max/num_k);  % (A^-1)
    disp(['-- k=' num2str(k(nk)) '/' num2str(k_max) ' A^-1']);

    % Transverse k (x,y)
    l = 1; m = 0;
    lm    = sqrt((l^2)+(m^2));
    kx    = (l/lm)*k(nk)*a*1e10;
    ky    = (m/lm)*k(nk)*a*1e10;
    kt    = (kx^2)+(ky^2);
    kz    = -i*(1/dz);

    for (zz = 1:length(z))
        disp(['z=' num2str(z(zz)) ' m']);
        
        % Luttinger Hamiltonian
        P  = Ev + C*g1_vec(zz)*(kt + kx^2);
        Q  = C*g2_vec(zz)*(kt - 2*kz^2);
        R  = -(sqrt(3)*C*g2_vec(zz)*((kx^2)-(ky^2)))+(i*2*C*g3_vec(zz)*sqrt(3)*kx*ky);
        S  = 2*C*g3_vec(zz)*sqrt(3)*((kx-(i*ky))*kz);

        V_z = V_vb(zz,2);
        
        % 4X4
        H4 = -[P+Q-V_z  -S       R       0;
                -S'    P-Q-V_z   0       R;
                 R'      0     P-Q-V_z   S;
                 0       R'      S'     P+Q-V_z];

        [V,D]  = eig(H4);
        eiglst = sum(D);
        
        F_E_v_4{zz,nk}.E = real(eiglst);
        F_E_v_4{zz,nk}.F = V;     
    end

    % 6X6
    %     H6  =  -[P+Q        -S            R           0       -S/sqrt(2)       sqrt(2)*R;
    %              -S'        P-Q           0           R       -sqrt(2)*Q       sqrt(1.5)*S;
    %               R'         0           P-Q          S        sqrt(1.5)*S'    sqrt(2)*Q;
    %               0          R'           S'         P+Q      -sqrt(2)*R'     -S'/sqrt(2);
    %           -S'/sqrt(2) -sqrt(2)*Q'  sqrt(1.5)*S -sqrt(2)*R     P+del             0;
    %            sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q'  -S/sqrt(2)       0             P+del];
    %
    %     [V,D]=eig(H6);
    %     eiglst = sum(D);
    %     ELK6(Nk,:) = sort(real(eiglst));

    % Get the physical state indices
    %[D_f, V_f] = SelectEigenvectors(D(1:Np),V(1:Np,1:Np));

%     D_f = D; V_f = V;
% 
%     E_v_all(nk,:)  = D_f(1:5);
%     wf_v_all{nk}   = NormaliseWfs(z,V_f(1:Np,1:5));
%     E_c_all(nk,:)  = ones(length(E_c),1).*Eg + E_c' + ...
%         ones(length(E_c),1).*(((hbar^2)*k2)./(2*m_eff_av*e_0*(a^2)));

    %     figure(5); hold on;
    %     for (mm=1:length(wf_v_all{nk}(:,1)))
    %         plot(abs(wf_v_all{nk}(mm,:)).^2);
    %     end
end

% % Grouping the band and sub-bands
% hh_subbands = [ 1, 1 ;
%     2, 3 ];
% lh_subbands = [ 1, 2 ;
%     2, 5 ];
% 
% %% Calculating the optical transmission strength
% %wf_c_new = wf_c/1e9;
% % for (nk=1:(num_k+1))
% %     wf_v = wf_v_all{nk}/1e10;
% %     for (index_vb=1:4)
% %
% %         for (index_vc=1:length(wf_c(:,1)))
% %             if (sum(index_vb==index_hh)>0)
% %                 overlap_hh = trapz(z, (wf_c_new(index_vc,:).*wf_v(index_vb,:)));
% %                 overlap_lh = 0;
% %             else
% %                 overlap_lh = trapz(z, (wf_c_new(index_vc,:).*wf_v(index_vb,:)));
% %                 overlap_hh = 0;
% %             end
% %
% %             rts_TE(index_vb,index_vc,nk) =  0.5*(abs(overlap_hh^2) + (1/3)*abs(overlap_lh^2));
% %             rts_TM(index_vb,index_vc,nk) =  (2/3)*abs(overlap_lh^2);
% %         end
% %
% %     end
% % end:
% 
% for (kk=1:num_k+1)
%     wf_v = wf_v_all{kk};
%     wf_v(isnan(wf_v)) = 0;
%     overlap_temp = zeros(length(wf_c(:,1)),4);
%     for (ii_cb=1:length(wf_c(:,1)))
%         for (ii_vb=1:4)
%             overlap_temp(ii_cb,ii_vb) = abs(trapz(z, (conj(wf_c(ii_cb,:)).*wf_v(ii_vb,:)))).^2;
%         end
%     end
%     overlap{kk} = overlap_temp;
% 
%     rts_TE(kk) = 0.5*(abs(overlap_temp(1)) + (1/3)*abs(overlap_temp(2)));
%     rts_TM(kk) = (2/3)*abs(overlap_temp(2));
% end
% 
% % for (kk=1:num_k+1)
% %     overlap_mat(kk,:) = overlap{kk};
% % end
% %
% 
% 
% %% Plotting
% 
% k     = k*10;     % (nm^-1)
% 
% % Potential profile
% figure(1); hold on;
% for (vb_index=1:length(E_v_all(1,:)))
%     plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, E_v_all(1,vb_index).*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), ':r');
%     text(z(Nb+Nw+10)/1e-10, E_v_all(1,vb_index), sprintf('h%d',vb_index));
% end
% 
% % Wavefunctions/probabity amplitude
% figure(2);
% wf_vb_0 = wf_v_all{1};
% subplot(312); hold on;
% for (ii=1:length(hh_subbands(:,2)))    % heavy holes
%     plot(z/1e-10, abs(wf_vb_0(hh_subbands(ii,2),:)).^2);
%     [max_wf_vb,max_index] = max(abs(wf_vb_0(hh_subbands(ii,2),:)).^2);
%     text(z(max_index)/1e-10, abs(wf_vb_0(hh_subbands(ii,2),max_index)).^2, sprintf('hh%d',hh_subbands(ii,1)));
% end
% ylabel('|F_h_h|^2'); box on; grid on;
% subplot(313); hold on;
% for (ii=1:length(lh_subbands(:,2)))    % light holes
%     plot(z/1e-10, abs(wf_vb_0(lh_subbands(ii,2),:)).^2);
%     [max_wf_vb,max_index] = max(abs(wf_vb_0(lh_subbands(ii,2),:)).^2);
%     text(z(max_index)/1e-10, abs(wf_vb_0(lh_subbands(ii,2),max_index)).^2, sprintf('lh%d',lh_subbands(ii,1)));
% end
% ylabel('|F_l_h|^2');
% grid on; box on;
% xlabel('z (Angstrom)');
% 
% % Dispersion
% figure(3);
% hold on;
% for (ii=1:length(E_v_all(1,:)))
%     plot(k, E_v_all(:,ii)*1e3);
% end
% for (jj=1:length(E_c_all(1,:)))
%     plot(k, E_c_all(:,jj)*1e3);
% end
% title(['GaAs/Ga_{' num2str(x) '}Al_{' num2str(1-x) '}As, W_G_a_A_s=' num2str(W) 'nm, Eg(k=0)=' num2str(Eg*1e3) 'meV: Conduction/Valence Bands E(k)']);
% xlabel(['k (nm^-^1) : <' num2str(l) num2str(m) '>']);
% ylabel('Energy (meV)')
% grid on;
% 
% %% Optical transitions
% 
% for (kk=1:num_k+1)
%     for (ii=1:length(overlap{1}(:,1)))
%         overlap_rts{ii}(kk,:) = overlap{kk}(ii,:);
%     end
% end
% 
% for (cc=1:length(overlap_rts))
%     figure; hold on; box on;
%     for(ii=1:length(overlap_rts{cc}(1,:)))
%         [max_over,index] = max(overlap_rts{cc}(:,ii));
%         plot(k, overlap_rts{cc});
%         legend(sprintf('cb%d-vb%d', cc, ii));
%         %text(k(index),max_over, ['\leftarrow' sprintf('cb%d-vb%d', cc, ii)] ,...
%         %     'HorizontalAlignment','left')
%     end
% end
