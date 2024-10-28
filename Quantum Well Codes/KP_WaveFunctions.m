
function [energies,wf,rts] = KP_WaveFunctions(x,W)

global Consts;
Constants;

%% 1. Luttinger-Kohn parameters

g1 = 6.85; g2 = 2.1; g3 = 2.9;    %GaAs

Gw1 = (Consts.hbar^2)*g1/(2*Consts.m_0);
Gw2 = (Consts.hbar^2)*g2/(2*Consts.m_0);
Gw3 = (Consts.hbar^2)*g3/(2*Consts.m_0);
Cw = Gw1+Gw2;
Dw = -2*Gw2+Gw1;
Aw = Gw1-Gw2;
Bw = 2*Gw2+Gw1;
Ew = Gw2+Gw3;

g1 = 3.45; g2 = 0.68; g3 = 1.29;  %AlAs

a1 = (Consts.hbar^2)*g1/(2*Consts.m_0);
Gb1 = ((1-x)*Gw1)+(x*a1);
a2 = (Consts.hbar^2)*g2/(2*Consts.m_0);
Gb2 = ((1-x)*Gw2)+(x*a2);
a3 = (Consts.hbar^2)*g3/(2*Consts.m_0);
Gb3 = ((1-x)*Gw3)+(x*a3);
Cb = Gb1+Gb2;
Db = -2*Gb2+Gb1;
Ab = Gb1-Gb2;
Bb = 2*Gb2+Gb1;
Eb = Gb2+Gb3;

Ev=0; Evb=(((1-x)*0)+(x*0.75))*Consts.e_0;
Eg = 1.426;   % eV


%% 2. Structure parameters

Nw = 40;
a = W/((Nw-1)*1e9);
%Nw = round(W/(a*1e9))+1
Nb = Nw;
Np = Nb+Nw+Nb;
Z=zeros(Np);

z = linspace(0,Np*a,Np);
dz = a;
dzz = dz^2;

A  = [Ab*ones(1,Nb),Aw*ones(1,Nw),Ab*ones(1,Nb)];
B  = [Bb*ones(1,Nb),Bw*ones(1,Nw),Bb*ones(1,Nb)];
C  = [Cb*ones(1,Nb),Cw*ones(1,Nw),Cb*ones(1,Nb)];
D  = [Db*ones(1,Nb),Dw*ones(1,Nw),Db*ones(1,Nb)];
E  = [Eb*ones(1,Nb),Ew*ones(1,Nw),Eb*ones(1,Nb)];
G2 = [Gb2*ones(1,Nb),Gw2*ones(1,Nw),Gb2*ones(1,Nb)];
G3 = [Gb3*ones(1,Nb),Gw3*ones(1,Nw),Gb3*ones(1,Nb)];
V  = [Evb*ones(1,Nb),zeros(1,Nw), Evb*ones(1,Nb)];

%% 3. Conduction band (using the shooting technique)

hs_structure = [Nb*a, Nb, x ; (W/1e9), Nw, 0 ; Nb*a, Nb, x];
[E_c,V_cb,V_vb,m_eff,con_profile] = OneParticleEnergies('e', hs_structure, 0, 2, Np, a);
wf_c  = OneParticleWavefunctions('e',E_c,V_cb,m_eff,con_profile);
wf_c  = wf_c./sqrt(trapz(z, wf_c(1,:).^2));
Evb = max(V_vb(:,2))./Consts.e_0;
V_vb = -V_vb(:,2)./Consts.e_0;

size(wf_c)


%% 4. Velence band (4X4 k.p)

for nk=1:100
    k(nk)=(nk-1)/500; % in A^-1
    l=0;m=1;lm=sqrt((l^2)+(m^2));
    kx=(l/lm)*k(nk)*1e10;ky=(m/lm)*k(nk)*1e10;
    k2=(kx^2)+(ky^2);
    k_t = sqrt(k2);
    k_t_vec(nk) = k_t;

    % Building the 2X2 diagonalized k.p Hamiltonian matrix ---------
    H = zeros(2*Np);
    for (zz = 1:Np)
        if (zz==1)
            H_diag  = [C(zz)*k_t^2+(D(zz+1)+2*D(zz))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                 ;
                                sqrt(3)*G2(zz)*k_t^2               , A(zz)*k_t^2+(B(zz+1)+2*B(zz))/(2*dzz)+V(zz)   ];

            H_off_p = [       -(D(zz+1)+D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ;
                        sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)        ];

            H_off_m = [      -(D(zz))/(2*dzz)        ,      sqrt(3)*k_t*(G3(zz))/(2*dz)       ;
                       -sqrt(3)*k_t*(G3(zz))/(2*dz)  ,            -(B(zz))/(2*dzz)            ];

            H(1:2,1:2) = H_diag;
            H(1:2,3:4) = H_off_p;
            H(3:4,1:2) = H_off_m;
        elseif (zz==Np)
            H_diag  = [C(zz)*k_t^2+(2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                       ;
                                sqrt(3)*G2(zz)*k_t^2               , A(zz)*k_t^2+(2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];

            H_off_p = [       -(D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz))/(2*dz)  ;
                        sqrt(3)*k_t*(G3(zz))/(2*dz)   ,            -(B(zz))/(2*dzz)        ];

            H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)        ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                       -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz) ,            -(B(zz)+B(zz-1))/(2*dzz)            ];

            H(2*Np-1:2*Np,2*Np-1:2*Np) = H_diag;
        else
            H_diag  = [C(zz)*k_t^2+(D(zz+1)+2*D(zz)+D(zz-1))/(2*dzz)+V(zz) ,          sqrt(3)*G2(zz)*k_t^2                       ;
                                sqrt(3)*G2(zz)*k_t^2                       , A(zz)*k_t^2+(B(zz+1)+2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ];

            H_off_p = [       -(D(zz+1)+D(zz))/(2*dzz)        ,     -sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ;
                        sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(2*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)        ];

            H_off_m = [      -(D(zz)+D(zz-1))/(2*dzz)         ,      sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ;
                       -sqrt(3)*k_t*(G3(zz)+G3(zz-1))/(2*dz)  ,            -(B(zz)+B(zz-1))/(2*dzz)            ];

            H(2*zz-1:2*zz,2*zz-1:2*zz)   = H_diag;
            H(2*zz-1:2*zz,2*zz+1:2*zz+2) = H_off_p;
            H(2*zz+1:2*zz+2,2*zz-1:2*zz) = H_off_m;
        end
    end

    % Calculate the eigenvectors of the Hamiltonian ----------------
    [nk sum(sum(abs(H-H')))];
    [F,En] = eig(H);
    En = diag(En)./Consts.e_0;
    [En,I] = sort(real(En));
    En = -En;
    F1{nk} = F(1:2:end,I);
    F2{nk} = F(2:2:end,I);

    E_v{nk} = En(1:4);

    E1(nk)=En(1); E2(nk)=En(2); E3(nk)=En(3); E4(nk)=En(4);

    for (ii=1:length(E_v{1}))
        F_U_norm = sqrt(trapz(z, F1{nk}(:,ii).^2+F2{nk}(:,ii).^2));
        Fhh{nk}(:,ii) = F1{nk}(:,ii)./F_U_norm;
        Flh{nk}(:,ii) = F2{nk}(:,ii)./F_U_norm;
    end

    % Overlap integrals --------------------------------------------
    for (cb_index=1:length(E_c))
        for (vb_index=1:length(E_v{1}))
            overlap1(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*Fhh{nk}(:,vb_index)).^2;
            overlap2(nk,cb_index,vb_index) = trapz(z, wf_c(cb_index,:)'.*Flh{nk}(:,vb_index)).^2;

            rtsTE(nk,cb_index,vb_index) = (3/4)*overlap1(nk, cb_index,vb_index) + ...
                                          (3/12)*overlap2(nk, cb_index,vb_index);
            rtsTE_k(vb_index,cb_index,nk) = rtsTE(nk,cb_index,vb_index);
            rtsTM(nk,cb_index,vb_index) = overlap2(nk, cb_index,vb_index);
            rtsTM_k(vb_index,cb_index,nk) = rtsTM(nk,cb_index,vb_index);
        end
    end
end

energies.E_c = E_c;
energies.E_v = E_v;
wf.wf_c = wf_c';
wf.wf_v.hh = Fhh;
wf.wf_v.lh = Flh;
rts.TE = rtsTE_k;
rts.TM = rtsTM_k;
wf.z = z;