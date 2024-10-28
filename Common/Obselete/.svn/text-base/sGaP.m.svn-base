% M file GaP.m
%*******************************************************
% M file with GaP data
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

EG = 2.26;					% Bandabstand bei 300 K in eV
epsrel = 11.2;				% relative DK
me = 0.58;					% effektive Masse des Leitungsbandes
mh = 0.54;					% effektive Masse des Valenzbandes
nue_e = 3;					% Zahl der Täler des Leitungsbandes
mye_i = 160;            % Elektronenbeweglichkeit im eigenleitenden Material bei 300 K
myh_i = 135;            % Löcherbeweglichkeit im eigenleitenden Material bei 300 K

ED = EH * me/epsrel^2;  % Bindungsenergie des Donators in eV
aBe = aB0 * epsrel/me;  % Störstellenradius des Donators in m
EA = EH * mh/epsrel^2;  % Bindungsenergie des Akzeptors in eV
aBh = aB0 * epsrel/mh;  % Störstellenradius des Akzeptors in m

% Berechnung für 300 K
Nc_0 = 2 * nue_e * me^1.5 * 1.254708e19;  % eff. Zustandsdichte des Leitungsbands in cm-3
Nv_0 = 2 * mh^1.5 * 1.254708e19;          % eff. Zustandsdichte des Valenzbands in cm-3
ni_square_0 = Nc_0 * Nv_0 * exp(-EG/kT); 
ni_0 = sqrt(ni_square_0);                  % intrinsische Ladungsträgerkonzentration in cm-3


