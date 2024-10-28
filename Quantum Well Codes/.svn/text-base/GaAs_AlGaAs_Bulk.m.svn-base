clear all; %close all;

z=zeros(5);Z=zeros(10);

%% Constants (all MKS, except energy which is in eV)

hbar=1.055e-34;q=1.602e-19;a=2.45e-10*4/sqrt(3);m=9.110e-31;

%% Valence band Luttinger-Kohn parameters (GaAs)

Ev   = 0;
del  = 0.3;

g1   = 6.85; g2 = 2.1; g3 = 2.9;
t1   = (hbar^2)*g1/(2*m*q*(a^2));
t2   = (hbar^2)*g2/(2*m*q*(a^2));
t3   = (hbar^2)*g3/(2*m*q*(a^2));

Nt   = 501;
kk   = 1*linspace(0,1,Nt);   % k*a
k_vec = [];
a     = 5.6605e-10;          % GaAs,AlAs lattice constant (meter)

l = 1;    m = 0;   n = 0;       %X-direction
%l = 0.5;  m = 0.5; n = 0.5;     %L-direction

%% Simulation

for Nk=1:Nt
    k = 2*pi*kk(Nk)*[l m n];
    k_vec = [k_vec, norm(k)];
    
    % Valence band Luttinger-Kohn model
    P  = Ev + (t1*sum(k.*k));
    Q  = t2*((k(1)^2) + (k(2)^2) - (2*(k(3)^2)));
    R  = -(sqrt(3)*t2*((k(1)^2)-(k(2)^2)))+(i*2*t3*sqrt(3)*k(1)*k(2));
    S  = 2*t3*sqrt(3)*((k(1)-(i*k(2)))*k(3));
    
    % 4X4 
    H4 = -[P+Q -S    R   0;
           -S' P-Q   0   R;
            R'  0   P-Q  S;
            0   R'   S' P+Q];
        
    [V,D]  = eig(H4);
    eiglst = sum(D);
    ELK4(Nk,:) = sort(real(eiglst));
    
    % 6X6
    H6  =  -[P+Q        -S            R           0       -S/sqrt(2)       sqrt(2)*R;
             -S'        P-Q           0           R       -sqrt(2)*Q       sqrt(1.5)*S;
              R'         0           P-Q          S        sqrt(1.5)*S'    sqrt(2)*Q;
              0          R'           S'         P+Q      -sqrt(2)*R'     -S'/sqrt(2);
          -S'/sqrt(2) -sqrt(2)*Q'  sqrt(1.5)*S -sqrt(2)*R     P+del             0;
           sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q'  -S/sqrt(2)       0             P+del];
       
    [V,D]=eig(H6);
    eiglst = sum(D);
    ELK6(Nk,:) = sort(real(eiglst));
end

%% Plotting

kk=-kk;        % L-direction - kk = -kk
k_vec = k_vec;

figure(1);
hold on;
plot(k_vec./a,ELK4,'b');
plot(k_vec./a,ELK6,'r');
xlabel('k (m^{-1}) ')
ylabel('Energy (eV) ')
grid on;

figure(2);
hold on;
plot((k_vec./a)/(2*pi/a),ELK4,'b');
plot((k_vec./a)/(2*pi/a),ELK6,'r');
xlabel('k/(2\pi/a)')
ylabel('Energy (eV) ')
grid on;
