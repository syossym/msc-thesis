% M file InP.m
%*******************************************************
% M file with InP data
%
% EG          bandgap at 300 K in eV
% epsrel      relative dielectric constant
% me          effektive mass of the conduction band
% mh          effektive mass of the valence band
% nue_e       number of conduction band minima
% ED          donor binding energy in eV
% aBe         donor electron radius in m
% EA          acceptor binding energy in eV
% aBh         acceptor hole radius in m
% Nc_0        effective density of conduction band states at 300 K in cm-3
% Nv_0        effective density of valence band states at 300 K in cm-3
% ni_0        intrinsic carrier concentration at 300 K in cm-3
% mye_i       electron mobility in intrinsic material at 300 K in cm^2/(Vs)
% myh_i       hole mobility in intrinsic material at 300 K in cm^2/(Vs)

konstanten;             % calls basic physical data

EG = 1.34;              % Bandabstand bei 300 K in eV
epsrel = 12.5;          % relative DK
me = 0.077;             % effektive Masse des Leitungsbandes
mh = 0.64;              % effektive Masse des Valenzbandes
nue_e = 1;					% Zahl der T�ler des Leitungsbandes
mye_i = 5370;           % Elektronenbeweglichkeit im eigenleitenden Material bei 300 K
myh_i = 150;            % L�cherbeweglichkeit im eigenleitenden Material bei 300 K

ED = EH * me/epsrel^2;  % Bindungsenergie des Donators in eV
aBe = aB0 * epsrel/me;  % St�rstellenradius des Donators in m
EA = EH * mh/epsrel^2;  % Bindungsenergie des Akzeptors in eV
aBh = aB0 * epsrel/mh;  % St�rstellenradius des Akzeptors in m

% Berechnung f�r 300 K
Nc_0 = 2 * nue_e * me^1.5 * 1.254708e19;  % eff. Zustandsdichte des Leitungsbands in cm-3
Nv_0 = 2 * mh^1.5 * 1.254708e19;          % eff. Zustandsdichte des Valenzbands in cm-3
ni_square_0 = Nc_0 * Nv_0 * exp(-EG/kT); 
ni_0 = sqrt(ni_square_0);                 % intrinsische Ladungstr�gerkonzentration in cm-3

