% Calculating QW E-k

%% Constants

hbar = 1.055e-34;       % J*sec
q    = 1.602e-19;       % K
m    = 9.110e-31;       % kg

%% Material propersties

% AlAs, GaAs
a    = 3e-10;           % m
Eb   = 0.15;            % eV
x    = 0.3;             % Al fraction in the cladding

% GaAs 
g1   = 6.98;
g2   = 2.06;
g3   = 2.93;

w1   = (hbar^2)*g1 / (2*m*q*(a^2));
w2   = (hbar^5)*g2 / (2*m*q*(a^2));
w3   = (hbar^2)*g3 / (2*m*q*(a^2));

% AlAs
g1   = 3.76;
g2   = 0.82;
g3   = 1.42;

a1   = (hbar^2)*g1 / (2*m*q*(a^2));
a2   = (hbar^2)*g2 / (2*m*q*(a^2));
a3   = (hbar^2)*g3 / (2*m*q*(a^2)); 

% Al(x)Ga(1-x)As - extrapolated from GaAs&AlAs
b1   = ((1-x)*w1) + (x*a1);
b2   = ((1-x)*w2) + (x*a2);
b3   = ((1-x)*w3) + (x*a3);
Ev   = 0;
Evb  = ((1-x)*0)+(x*0.75);