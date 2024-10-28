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
w1  = (hbar^2)*g1/(2*m_0*e_0*(a^2));
w2  = (hbar^2)*g2/(2*m_0*e_0*(a^2));
w3  = (hbar^2)*g3/(2*m_0*e_0*(a^2));
Ev  = 0;
Eg  = 1.426;                % band gap (eV)

% AlAs
g1  = 3.76; g2 = 0.82; g3 = 1.42;
a1  = (hbar^2)*g1/(2*m_0*e_0*(a^2));
a2  = (hbar^2)*g2/(2*m_0*e_0*(a^2));
a3  = (hbar^2)*g3/(2*m_0*e_0*(a^2));

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
W        = 9; %num_mono*W_mono/10;   % well width (nm)
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

%hs_structure = [10*(Nb-1)*(a*1e9), Nb, x ; 10*W, Nw, 0 ; 10*(Nb-1)*(a*1e9), Nb, x];   % the width are given in A
hs_structure                        = [10*Wb, Nb, x ; 10*W, Nw, 0 ; 10*Wb, Nb, x];
[E_c,V_cb,V_vb,m_eff,con_profile]   = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
wf_c                                = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
m_eff_av                            = mean(m_eff(:,2));
z                                   = V_cb(:,1);
Evb = max(V_vb(:,2))/e_0;
wf_c = NormaliseWfs(z, wf_c.');

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
num_k = 10;     % size of k-axis

for (nk=1:(num_k+1))

    k(nk) = (nk-1)*(k_max/num_k);  % (A^-1)
    disp(['k=' num2str(k(nk)) '/' num2str(k_max) ' A^-1']);

    % Transverse k (x,y)
    l = 1; m = 1;
    lm  = sqrt((l^2)+(m^2));
    kx  = (l/lm)*k(nk)*a*1e10;
    ky  = (m/lm)*k(nk)*a*1e10;
    k2  = (kx^2)+(ky^2);
    k_t(nk) = sqrt(k2);

    % Conversion matrix
    if (nk==1)
        fi = 0;
    else
        fi    = atan(ky/kx);
    end
    alpha = (1/sqrt(2))*exp(i*(3*pi/4-3*fi/2));
    beta  = (1/sqrt(2))*exp(-i*(pi/4-fi/2));
    U_t   = [conj(alpha)   0      0      -alpha ;
                 0      -beta  conj(beta)   0   ;
                 0       beta  conj(beta)   0   ;
              conj(alpha)   0     0        alpha ];

    % Building the 4X4 Luttinger Hamiltonian for the entire structure
    t     = [b1*ones(1,Nb) w1*ones(1,Nw-1) b1*ones(1,Nb)];
    tt    = [0 t]+[t 0];
    Ebk   = Evb+(b1*k2);
    Ewk   = (w1*k2);
    Ebwk  = (Ebk+Ewk)/2;
    U     = Ev+[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    P     = -diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    t     = -2*[b2*ones(1,Nb) w2*ones(1,Nw-1) b2*ones(1,Nb)];
    tt    = [0 t]+[t 0];
    Ebk   = b2*k2;
    Ewk   = w2*k2;
    Ebwk  = (Ebk+Ewk)/2;
    U     = [Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    Q     = -diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    Ebk   = -(sqrt(3)*b2*((kx^2)-(ky^2)))+(i*2*b3*sqrt(3)*kx*ky);
    Ewk   = -(sqrt(3)*w2*((kx^2)-(ky^2)))+(i*2*w3*sqrt(3)*kx*ky);
    Ebwk  = (Ebk+Ewk)/2;
    U     = [Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    R     = diag(U);

    t     = -2*i*sqrt(3)*(kx-(i*ky))*[b3*ones(1,Nb) w3*ones(1,Nw-1) b3*ones(1,Nb)]/2;
    S     = diag(t,1)-diag(t,-1);

    H     = [P+Q Z;
             Z P+Q];

    HL    = [P-Q Z;
             Z P-Q];

    HC    = [-S  R;
              R' S'];

    H     = -[H   HC;
              HC' HL];

        
    % The energies and states
    [V,D] = eig(H);
    D     = diag(D);
    [D,I] = sort(real(-D));
    D     = -(D)';
    V     = V(:,I);

    % Envelope functions
    F1{nk} = V(1:Np, 1:10);
    F2{nk} = V(Np+1:2*Np, 1:10);
    F3{nk} = V(2*Np+1:3*Np, 1:10);
    F4{nk} = V(3*Np+1:4*Np, 1:10);
    F_norm = sqrt(trapz(z,F1{nk}.^2+F2{nk}.^2+F3{nk}.^2+F4{nk}.^2,1));
    
    for(jj = 1:length(F1{nk}(1,:)))

        % Transform the vectors
        F       = [F1{nk}(:,jj).'; F2{nk}(:,jj).'; F3{nk}(:,jj).'; F4{nk}(:,jj).'];
        for (ff=1:Np)
           F_trans(:,ff) = U_t*F(:,ff); 
        end
        
        %F_trans = U_t*F;

        F1{nk}(:,jj) = F1{nk}(:,jj)./F_norm(jj);
        F2{nk}(:,jj) = F2{nk}(:,jj)./F_norm(jj);
        F3{nk}(:,jj) = F3{nk}(:,jj)./F_norm(jj);
        F4{nk}(:,jj) = F4{nk}(:,jj)./F_norm(jj);

        F_trans_norm_U = sqrt(trapz(z,F_trans(1,:).^2+F_trans(2,:).^2));
        F_trans_norm_D = sqrt(trapz(z,F_trans(3,:).^2+F_trans(4,:).^2));
        F1_u{nk}(:,jj) = F_trans(1,:).'./F_trans_norm_U;
        F2_u{nk}(:,jj) = F_trans(2,:).'./F_trans_norm_U;
        F3_d{nk}(:,jj) = F_trans(3,:).'./F_trans_norm_D;
        F4_d{nk}(:,jj) = F_trans(4,:).'./F_trans_norm_D;
        
        % Overlap integrals
        overlap_c1_F1_u{nk}(:,jj) = trapz(z, abs(wf_c(1,:).').*abs(F1_u{nk}(:,jj)));
        overlap_c2_F1_u{nk}(:,jj) = trapz(z, abs(wf_c(2,:).').*abs(F1_u{nk}(:,jj)));
        overlap_c1_F2_u{nk}(:,jj) = trapz(z, abs(wf_c(1,:).').*abs(F2_u{nk}(:,jj)));
        overlap_c2_F2_u{nk}(:,jj) = trapz(z, abs(wf_c(2,:).').*abs(F2_u{nk}(:,jj)));
    end
    
    % Get the physical state indices
    %[D_f, V_f] = SelectEigenvectors(D(1:Np),V(1:Np,1:Np));

    D_f = D; V_f = V;

    E_v_all(nk,:)  = 10*D_f(1:10);
    %wf_v_all{nk}   = NormaliseWfs(z,V_f(1:Np,1:5));
    E_c_all(nk,:)  = ones(length(E_c),1).*Eg + E_c' + ...
        ones(length(E_c),1).*(((hbar^2)*k2)./(2*m_eff_av*e_0*(a^2)));

    % Plotting
    figure(2);

    subplot(312); box on; grid on;
    for (cb_index=1:length(E_c))        % electrons
        hold on; 
        plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
        [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
        text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
        hold off;
    end
    for (ii=1:length(F1_u{nk}(1,:)))    % heavy holes
        hold on; 
        plot(z/1e-10, abs(F1_u{nk}(:,ii)).^2);
        [m,max_index] = max(abs(F1_u{nk}(:,ii)).^2);
        text(z(max_index)/1e-10, abs(F1_u{nk}(max_index,ii)).^2, sprintf('h%d',ii));
        hold off;
    end
    ylabel('|F^U_1|^2');

%     subplot(313); box on; grid on;
%     for (cb_index=1:length(E_c))        % electrons
%         hold on;
%         plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
%         [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
%         text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
%         hold off;
%     end
%     for (ii=1:length(F2_u{nk}(1,:)))    % light holes
%         hold on;
%         plot(z/1e-10, abs(F2_u{nk}(:,ii)).^2);
%         [m,max_index] = max(abs(F2_u{nk}(:,ii)).^2);
%         text(z(max_index)/1e-10, abs(F2_u{nk}(max_index,ii)).^2, sprintf('h%d',ii));
%         hold off;
%     end
%     ylabel('|F^U_2|^2');
%     xlabel('z (Angstrom)');
%     
%     figure(5);
%     subplot(411); box on; grid on;
%     for (cb_index=1:length(E_c))        % electrons
%         hold on;
%         plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
%         [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
%         text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
%         hold off;
%     end
%     for (ii=1:10)
%         hold on;
%         plot(z/1e-10, abs(F1{nk}(:,ii)).^2);
%     end
%     ylabel('|F_1|^2');
%     subplot(412); box on; grid on;
%     for (cb_index=1:length(E_c))        % electrons
%         hold on;
%         plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
%         [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
%         text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
%         hold off;
%     end
%     for (ii=1:10)
%         hold on;
%         plot(z/1e-10, abs(F2{nk}(:,ii)).^2);
%     end
%     ylabel('|F_2|^2');
%     subplot(413); box on; grid on;
%     for (cb_index=1:length(E_c))        % electrons
%         hold on;
%         plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
%         [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
%         text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
%         hold off;
%     end
%     for (ii=1:10)
%         hold on;
%         plot(z/1e-10, abs(F3{nk}(:,ii)).^2);
%     end
%     ylabel('|F_3|^2');
%     subplot(414);  box on; grid on;
%     for (cb_index=1:length(E_c))        % electrons
%         hold on;
%         plot(z/1e-10, abs(wf_c(cb_index,:)).^2, 'r');
%         [max_wf_c,max_index] = max(abs(wf_c(cb_index,:)).^2);
%         text(z(max_index)/1e-10, abs(wf_c(cb_index,max_index)).^2, sprintf('e%d',cb_index));
%         hold off;
%     end
%     for (ii=1:10)
%         hold on;
%         plot(z/1e-10, abs(F4{nk}(:,ii)).^2);
%     end
%     ylabel('|F_4|^2');
%     xlabel('z (Angstrom)');
end

clear H V D;

% Grouping the band and sub-bands
hh_subbands = [ 1, 1 ;
                2, 3 ];
lh_subbands = [ 1, 2 ;
                2, 5 ];

%wf_c_new = wf_c/1e9;
% for (nk=1:(num_k+1))
%     wf_v = wf_v_all{nk}/1e10;
%     for (index_vb=1:4)
%
%         for (index_vc=1:length(wf_c(:,1)))
%             if (sum(index_vb==index_hh)>0)
%                 overlap_hh = trapz(z, (wf_c_new(index_vc,:).*wf_v(index_vb,:)));
%                 overlap_lh = 0;
%             else
%                 overlap_lh = trapz(z, (wf_c_new(index_vc,:).*wf_v(index_vb,:)));
%                 overlap_hh = 0;
%             end
%
%             rts_TE(index_vb,index_vc,nk) =  0.5*(abs(overlap_hh^2) + (1/3)*abs(overlap_lh^2));
%             rts_TM(index_vb,index_vc,nk) =  (2/3)*abs(overlap_lh^2);
%         end
%
%     end
% end:

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

% for (kk=1:num_k+1)
%     overlap_mat(kk,:) = overlap{kk};
% end
%


%% Plotting

k     = k*10;     % (nm^-1)

% Potential profile
figure(1); hold on;
for (vb_index=1:length(E_v_all(1,:)))
    plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, E_v_all(1,vb_index).*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), ':r');
    text(z(Nb+Nw+10)/1e-10, E_v_all(1,vb_index), sprintf('h%d',vb_index));
end

% Wavefunctions/probabity amplitude
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
ylabel('Energy (meV)')
grid on;

%% Optical transitions

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

overlap = cell(2);

for (nk=1:(num_k+1))
   overlap{1,1} = [overlap{1,1}; overlap_c1_F1_u{nk}];
   overlap{1,2} = [overlap{1,2}; overlap_c1_F2_u{nk}];
   overlap{2,1} = [overlap{2,1}; overlap_c2_F1_u{nk}];
   overlap{2,2} = [overlap{2,2}; overlap_c2_F2_u{nk}];
end

for (cb_index=1:2)
    rts_TE{cb_index} = 0.5.*abs(overlap{cb_index,1}).^2 + (1/3)*abs(overlap{cb_index,2}).^2;
    rts_TM{cb_index} = (1/3)*abs(overlap{cb_index,2}).^2;

    
    figure(4);
    for (ii=1:4)  
        subplot(211);hold on;
        plot(k_t/(2*pi/(a*1e10)), rts_TE{cb_index}(:,ii));
        subplot(212);hold on;
        plot(k_t/(2*pi/(a*1e10)), rts_TM{cb_index}(:,ii));
    end

end

