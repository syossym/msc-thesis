
clear all; close all; clc;

%Constants (all MKS, except energy which is in eV)
hbar=1.055e-34;q=1.602e-19;a=3e-10;m=9.110e-31;

%Luttinger-Kohn parameters
g1=6.85;g2=2.1;g3=2.9;%GaAs
Gw1=(hbar^2)*g1/(2*m);
Gw2=(hbar^2)*g2/(2*m);
Gw3=(hbar^2)*g3/(2*m);
Aw = Gw1+Gw2; Bw = 2*Gw2-Gw1; Cw = Gw1-Gw2; Dw = 2*Gw2+Gw1; Ew = Gw2+Gw3;

g1=3.45;g2=0.68;g3=1.29;%AlAs
a1=(hbar^2)*g1/(2*m);
Gb1=(.7*Gw1)+(.3*a1);
a2=(hbar^2)*g2/(2*m);
Gb2=(.7*Gw2)+(.3*a2);
a3=(hbar^2)*g3/(2*m);
Gb3=(.7*Gw3)+(.3*a3);
Ab = Gb1+Gb2; Bb = 2*Gb2-Gb1; Cb = Gb1-Gb2; Db = 2*Gb2+Gb1; Eb = Gb2+Gb3;

Ev=0;Evb=((0.7*0)+(0.3*0.75))*q; 

Nw=18;Nb=Nw;Np=Nb+Nw+Nb;W=(Nw-1)*a*1e9,Z=zeros(Np);
dz = a;
dzz = dz^2;

A  = [Ab*ones(1,Nb),Aw*ones(1,Nw),Ab*ones(1,Nb)];
B  = [Bb*ones(1,Nb),Bw*ones(1,Nw),Bb*ones(1,Nb)];
C  = [Cb*ones(1,Nb),Cw*ones(1,Nw),Cb*ones(1,Nb)];
D  = [Db*ones(1,Nb),Dw*ones(1,Nw),Db*ones(1,Nb)];
E  = [Eb*ones(1,Nb),Ew*ones(1,Nw),Eb*ones(1,Nb)];
G3 = [Gb3*ones(1,Nb),Gw3*ones(1,Nw),Gb3*ones(1,Nb)];
V  = [Evb*ones(1,Nb),zeros(1,Nw), Evb*ones(1,Nb)];

% A  = [Ab*ones(1,Nb),0.5*(Ab+Aw),Aw*ones(1,Nw-2),0.5*(Ab+Aw),Ab*ones(1,Nb)];
% B  = [Bb*ones(1,Nb),0.5*(Bb+Bw),Bw*ones(1,Nw-2),0.5*(Bb+Bw),Bb*ones(1,Nb)];
% C  = [Cb*ones(1,Nb),0.5*(Cb+Cw),Cw*ones(1,Nw-2),0.5*(Cb+Cw),Cb*ones(1,Nb)];
% D  = [Db*ones(1,Nb),0.5*(Db+Dw),Dw*ones(1,Nw-2),0.5*(Db+Dw),Db*ones(1,Nb)];
% E  = [Eb*ones(1,Nb),0.5*(Eb+Ew),Ew*ones(1,Nw-2),0.5*(Eb+Ew),Eb*ones(1,Nb)];
% G3 = [Gb3*ones(1,Nb),0.5*(Gb3+Gw3),Gw3*ones(1,Nw-2),0.5*(Gb3+Gw3),Gb3*ones(1,Nb)];
% V  = [Evb*ones(1,Nb), Evb/2, zeros(1,Nw-2), Evb/2, Evb*ones(1,Nb)];

for nk=1:100
    k(nk)=(nk-1)/500;% in A^-1
    l=0;m=1;lm=sqrt((l^2)+(m^2));
    kx=(l/lm)*k(nk)*1e10;ky=(m/lm)*k(nk)*1e10;
    k2=(kx^2)+(ky^2);
    k_t = sqrt(k2)
    k_t_vec(nk) = k_t;

    % Building the matrix
    H = zeros(2*Np);
    for (zz = 1:Np)
        
        if (zz==1)
            H_diag  = [A(zz)*k_t^2-(B(zz+1)+2*B(zz))/(2*dzz)+V(zz) ,          0.5*sqrt(3)*E(zz)*k_t^2             ;
                                  0.5*sqrt(3)*E(zz)*k_t^2          , A(zz)*k_t^2+(B(zz+1)+2*B(zz))/(2*dzz)+V(zz)   ];
                      
            H_off_p = [       (B(zz+1)+B(zz))/(2*dzz)            ,     -0.5*sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(4*dz)  ;
                       0.5*sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(4*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)             ];
                             
            H_off_m = [       (B(zz))/(2*dzz)            ,      0.5*sqrt(3)*k_t*(G3(zz))/(4*dz)  ;
                       -0.5*sqrt(3)*k_t*(G3(zz))/(4*dz)  ,            -(B(zz))/(2*dzz)            ]; 
                   
            H(1:2,1:2) = H_diag;
            H(1:2,3:4) = H_off_p;
            H(3:4,1:2) = H_off_m;
        elseif (zz==Np)
            H_diag  = [A(zz)*k_t^2-(B(zz-1)+2*B(zz))/(2*dzz)+V(zz) ,          0.5*sqrt(3)*E(zz)*k_t^2             ;
                                  0.5*sqrt(3)*E(zz)*k_t^2          , A(zz)*k_t^2+(B(zz-1)+2*B(zz))/(2*dzz)+V(zz)   ];  
                             
            H_off_p = [       (B(zz))/(2*dzz)           ,     -0.5*sqrt(3)*k_t*(G3(zz))/(4*dz)  ;
                       0.5*sqrt(3)*k_t*(G3(zz))/(4*dz)  ,            -(B(zz))/(2*dzz)            ];
                   
            H_off_m = [       (B(zz-1)+B(zz))/(2*dzz)            ,     0.5*sqrt(3)*k_t*(G3(zz-1)+G3(zz))/(4*dz)  ;
                       -0.5*sqrt(3)*k_t*(G3(zz-1)+G3(zz))/(4*dz) ,            -(B(zz-1)+B(zz))/(2*dzz)             ];
                   
            H(2*Np-1:2*Np,2*Np-1:2*Np) = H_diag;
        else
            H_diag  = [A(zz)*k_t^2-(B(zz+1)+2*B(zz)+B(zz-1))/(2*dzz)+V(zz) ,                  0.5*sqrt(3)*E(zz)*k_t^2             ;
                                  0.5*sqrt(3)*E(zz)*k_t^2                  , A(zz)*k_t^2+(B(zz+1)+2*B(zz)+B(zz-1))/(2*dzz)+V(zz)   ];
                             
            H_off_p = [       (B(zz+1)+B(zz))/(2*dzz)            ,     -0.5*sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(4*dz)  ;
                       0.5*sqrt(3)*k_t*(G3(zz+1)+G3(zz))/(4*dz)  ,            -(B(zz+1)+B(zz))/(2*dzz)             ];
                   
            H_off_m = [       (B(zz-1)+B(zz))/(2*dzz)            ,     0.5*sqrt(3)*k_t*(G3(zz-1)+G3(zz))/(4*dz)  ;
                       -0.5*sqrt(3)*k_t*(G3(zz-1)+G3(zz))/(4*dz) ,            -(B(zz-1)+B(zz))/(2*dzz)             ];
                   
            H(2*zz-1:2*zz,2*zz-1:2*zz) = H_diag;
            H(2*zz-1:2*zz,2*zz+1:2*zz+2) = H_off_p;
            H(2*zz+1:2*zz+2,2*zz-1:2*zz) = H_off_m;
        end
        
        
    end
    [nk sum(sum(abs(H-H')))]
    [F,En]=eig(-H);
    En=diag(En)./q;
    En=(sort(real(En)))';
    E1(nk)=En(1);E2(nk)=En(2);E3(nk)=En(3);E4(nk)=En(4);
    E5(nk)=En(5);E6(nk)=En(6);E7(nk)=En(7);E8(nk)=En(8);
    E9(nk)=En(9);E10(nk)=En(10);E11(nk)=En(11);E12(nk)=En(12);
       
end

figure(2);
k=k*10;%per Angstrom to per nm
hold on
%h=plot(W,Ean1,'b');
%h=plot(W,Ean2,'b');
h=plot(k,E1,'b');
set(h,'linewidth',[2.0])
h=plot(k,E2,'bx');
h=plot(k,E3,'b');
h=plot(k,E4,'b+');
h=plot(k,E5,'b');
set(h,'linewidth',[2.0])
h=plot(k,E6,'bx');
h=plot(k,E7,'b');
h=plot(k,E8,'b+');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel(' k ( /nm ) ---> ')
ylabel(' Energy ( eV ) ---> ')
%axis([0 .5 -.1 0])
grid on