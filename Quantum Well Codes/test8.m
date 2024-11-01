
clear all; close all; clc; warning off;

%Constants (all MKS, except energy which is in eV)
hbar=1.055e-34;q=1.602e-19;a=3e-10;m=9.110e-31;

% Plotting definitions
h_struct = figure('Name','Structure');
h_env_func = figure('Name','Envelope Functions');
h_disp = figure('Name', 'Dispersion');
h_rts = figure('Name', 'Transition Matrix');

%Luttinger-Kohn parameters
g1=6.85;g2=2.1;g3=2.9;%GaAs
w1=(hbar^2)*g1/(2*m*q*(a^2));
w2=(hbar^2)*g2/(2*m*q*(a^2));
w3=(hbar^2)*g3/(2*m*q*(a^2));
g1=3.45;g2=0.68;g3=1.29;%AlAs
a1=(hbar^2)*g1/(2*m*q*(a^2));b1=(.7*w1)+(.3*a1);
a2=(hbar^2)*g2/(2*m*q*(a^2));b2=(.7*w2)+(.3*a2);
a3=(hbar^2)*g3/(2*m*q*(a^2));b3=(.7*w3)+(.3*a3);
Eg = 1.426; Ev=0;Evb=(0.7*0)+(0.3*0.75)

Nw=35;Nb=Nw;Np=Nb+Nw+Nb;W=(Nw-1)*a*1e9,Z=zeros(Np);
z = linspace(0,Np*a, Np);

hs_structure  = [Nb*a, Nb, 0.3 ; (W/1e9), Nw, 0 ; Nb*a, Nb, 0.3];
[E_c,V_cb,V_vb,m_eff,con_profile] = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
wf_c   = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
wf_c   = wf_c./sqrt(trapz(z, wf_c(1,:).^2)); %NormaliseWfs(z, wf_c.');
%V_vb = [-Evb*ones(1,Nb), -Evb/2, zeros(1,Nw-2), -Evb/2, -Evb*ones(1,Nb)];
Evb = max(V_vb(:,2))./q;
V_vb = -V_vb(:,2)./q;

% Plot the structure
figure(h_struct);
plot(z/1e-10, Eg*ones(1,Np) + V_cb(:,2)'./q, 'g', z/1e-10, V_vb, 'g', 'linewidth', 2);

for nk=1:100
    k(nk)=(nk-1)/1000;% in A^-1
    l=0;m=1;lm=sqrt((l^2)+(m^2));
    kx=(l/lm)*k(nk)*a*1e10;ky=(m/lm)*k(nk)*a*1e10;
    k2=(kx^2)+(ky^2);
    k_t = sqrt(k2);
    k_t_vec(nk) = k_t;

    t=[b1*ones(1,Nb) w1*ones(1,Nw-1) b1*ones(1,Nb)];tt=[0 t]+[t 0];
    Ebk=Evb+(b1*k2);Ewk=(w1*k2);Ebwk=(Ebk+Ewk)/2;
    U=Ev+[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    P=-diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    t=-2*[b2*ones(1,Nb) w2*ones(1,Nw-1) b2*ones(1,Nb)];tt=[0 t]+[t 0];
    Ebk=b2*k2;Ewk=w2*k2;Ebwk=(Ebk+Ewk)/2;
    U=[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    Q=-diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    Ebk=-(sqrt(3)*b2*((kx^2)-(ky^2)))+(i*2*b3*sqrt(3)*kx*ky);
    Ewk=-(sqrt(3)*w2*((kx^2)-(ky^2)))+(i*2*w3*sqrt(3)*kx*ky);
    Ebwk=(Ebk+Ewk)/2;
    U=[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    R=diag(U);

    t=-2*i*sqrt(3)*(kx-(i*ky))*[b3*ones(1,Nb) w3*ones(1,Nw-1) b3*ones(1,Nb)]/2;
    S=diag(t,1)-diag(t,-1);

    H=[P+Q Z;Z P+Q];HL=[P-Q Z;Z P-Q];
    HC=[-S R;R' S'];
    H=-[H HC;HC' HL];
    %[nk sum(sum(abs(H-H')))]

    [V,D]=eig(H);D=diag(D);
    [D,I]=sort(real(-D)); 
    D = -D';
    V = V(:,I);
    E1(nk)=D(1);E2(nk)=D(2);E3(nk)=D(3);E4(nk)=D(4);
    E5(nk)=D(5);E6(nk)=D(6);E7(nk)=D(7);E8(nk)=D(8);
    E9(nk)=D(9);E10(nk)=D(10);E11(nk)=D(11);E12(nk)=D(12);
    
    E_v{nk} = D(2:2:8);
    
    F1{nk} = V(1:Np,2:2:8);
    F2{nk} = V(Np:2*Np-1,2:2:8);
    F3{nk} = V(2*Np:3*Np-1,2:2:8);
    F4{nk} = V(3*Np:4*Np-1,2:2:8);
    
    for (ii=1:4)
        F_norm = sqrt(trapz(z, F1{nk}(:,ii).^2+F2{nk}(:,ii).^2+F3{nk}(:,ii).^2+F4{nk}(:,ii).^2));

        F1{nk}(:,ii) = F1{nk}(:,ii)./F_norm;
        F2{nk}(:,ii) = F2{nk}(:,ii)./F_norm;
        F3{nk}(:,ii) = F3{nk}(:,ii)./F_norm;
        F4{nk}(:,ii) = F4{nk}(:,ii)./F_norm;
    end
    
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
          
    for (ii=1:4)
       F_trans = U_t*[F1{nk}(:,ii).'; F2{nk}(:,ii).'; F3{nk}(:,ii).'; F4{nk}(:,ii).'];
       F_U_norm = sqrt(trapz(z, F_trans(1,:).^2+F_trans(2,:).^2));
       F1_U{nk}(:,ii) = -F_trans(1,:).'./F_U_norm;%./sqrt(trapz(z, F_trans(1,:).^2))
       F2_U{nk}(:,ii) = -F_trans(2,:).'./F_U_norm;%./sqrt(trapz(z, F_trans(2,:).^2));
    end
    
    temp = F2_U{nk}(:,2);
    F2_U{nk}(:,2) = F1_U{nk}(:,2);
    F1_U{nk}(:,2) = temp;
    
    temp = F1_U{nk}(:,4);
    F1_U{nk}(:,4) = -F2_U{nk}(:,4);
    F2_U{nk}(:,4) = temp;
   
    % Calculate the overlap integrals
    for (cb_index=1:length(E_c))
       for (vb_index=1:length(E_v{1}))
          overlap1(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*F1_U{nk}(:,vb_index)).^2;
          overlap2(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*F2_U{nk}(:,vb_index)).^2;
       
          rtsTE(nk,cb_index,vb_index) = (3/4)*overlap1(nk, cb_index,vb_index) + ... 
                                         (3/12)*overlap2(nk, cb_index,vb_index);
          rtsTM(nk,cb_index,vb_index) = overlap2(nk, cb_index,vb_index);
       end
    end
    
    % Plot the envelope functions
    figure(h_env_func);
    clf;
    colors = [1 0 0; 1 1 0; 0 1 0; 0 1 1];
    for(ii=1:4)
        subplot(211); grid on; box on;
        hold on; h_F1(ii) = plot(z/1e-10, F1_U{nk}(:,ii), 'Color', colors(ii,:)); hold off; 
        l_F1{ii} = sprintf('F_1^U - h%d',ii);  ylabel('F_1^u');
        subplot(212); grid on; box on;
        hold on; h_F2(ii) = plot(z/1e-10, F2_U{nk}(:,ii), 'Color', colors(ii,:)); hold off; 
        l_F2{ii} = sprintf('F_2^U - h%d',ii); ylabel('F_2^u');     
    end
    subplot(211); 
    hold on; plot(z/1e-10, wf_c, 'r--', 'linewidth', 2); hold off;
    legend(h_F1,l_F1); xlabel('z (A)'); 
    title(['Envelope function amplitudes, k_t=',num2str(k_t),' A^-^1']);
    subplot(212); 
    hold on; plot(z/1e-10, wf_c, 'r--', 'linewidth', 2); hold off;
    legend(h_F2,l_F2); xlabel('z (A)'); 
    %pause;
end

% Plot the k_t=0 bands and envelope functions
figure(h_struct); 
for(cb_index=1:length(E_c))
    hold on;
    plot(z(V_cb(:,2)==min(V_cb(:,2)))/1e-10, Eg + E_c(cb_index)*ones(1,length(z(V_cb(:,2)==min(V_cb(:,2))))), 'r', 'linewidth', 2 );
    text(z(Nb+Nw+10)/1e-10, Eg + E_c(cb_index), sprintf('e%d',cb_index));
    plot(z/1e-10, (Eg + E_c(cb_index))*ones(1,Np) + wf_c(cb_index,:)./1e5, 'r');
end
for(vb_index=1:length(E_v{1}))
    hold on;
    plot(z(V_vb==max(V_vb))/1e-10, E_v{1}(vb_index)*ones(1,length(z(V_vb==max(V_vb)))), 'b', 'linewidth', 2);
    text(z(Nb+Nw+10)/1e-10, E_v{1}(vb_index), sprintf('h%d',vb_index));
    plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np) + F1_U{1}(:,vb_index).'./1e5, ':b');
    plot(z/1e-10, E_v{1}(vb_index)*ones(1,Np) + F2_U{1}(:,vb_index).'./1e5, ':r');
end

% Plot the band dispersion relation
figure(h_disp);
k=k*10;%per Angstrom to per nm
hold on
%h=plot(W,Ean1,'b');
%h=plot(W,Ean2,'b');
h=plot(k_t_vec,E1,'b');
%set(h,'linewidth',[2.0])
h=plot(k_t_vec,E2,'bx');
h=plot(k_t_vec,E3,'b');
h=plot(k_t_vec,E4,'b+');
h=plot(k_t_vec,E5,'b');
%set(h,'linewidth',[2.0])
h=plot(k_t_vec,E6,'bx');
h=plot(k_t_vec,E7,'b');
h=plot(k_t_vec,E8,'b+');
%set(h,'linewidth',[2.0])
%set(gca,'Fontsize',[24])
xlabel('k_t (A^-^1)')
ylabel('Energy (eV) ')
%axis([0 .5 -.1 0])
grid on;

figure(h_rts);
colors = [1 0 0; 1 1 0; 0 1 0; 0 1 1];
for (jj=1:length(E_v{1}))  
    l_tran1{jj} = sprintf('c1-h%d',jj); 
    l_tran2{jj} = sprintf('c2-h%d',jj); 
    
    subplot(221); hold on; grid on; box on;
    h_tran1(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTE(:,1,jj)),30,'moving'), 'Color', colors(jj,:));
    %text(k_t_vec(end)*1e10, rtsTE(end,1,jj), sprintf('e1-h%d',jj));  
    subplot(222); hold on; grid on; box on;
    h_tran2(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTM(:,1,jj)),30,'moving'), 'Color', colors(jj,:));
    %text(k_t_vec(end)*1e10, rtsTM(end,1,jj), sprintf('e1-h%d',jj));
    subplot(223); hold on; grid on; box on;
    h_tran3(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTE(:,2,jj)),30,'moving'), 'Color', colors(jj,:));
    %text(k_t_vec(end)*1e10, rtsTE(end,2,jj), sprintf('e2-h%d',jj));
    subplot(224); hold on; grid on; box on;
    h_tran4(jj) = plot(k_t_vec*1e10, smooth(squeeze(rtsTM(:,2,jj)),30,'moving'), 'Color', colors(jj,:));
    %plot(k_t_vec*1e10, spline(k_t_vec*1e10, squeeze(rtsTM(:,2,jj))));
    %text(k_t_vec(end)*1e10, rtsTM(end,2,jj), sprintf('e2-h%d',jj));
end
subplot(221); legend(h_tran1, l_tran1); title('TE'); xlabel('k_t (m^-^1)')
subplot(222); legend(h_tran2, l_tran1); title('TM'); xlabel('k_t (m^-^1)')
subplot(223); legend(h_tran3, l_tran2); title('TE'); xlabel('k_t (m^-^1)')
subplot(224); legend(h_tran4, l_tran2); title('TM'); xlabel('k_t (m^-^1)');