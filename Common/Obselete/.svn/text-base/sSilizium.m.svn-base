% M file Silizium.m
%*******************************************************
% M file with silicon data
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

EG = 1.12;              % bandgap at 300 K in eV
epsrel = 11.4;				% relative dielectric constant
me = 0.32;					% effektive mass of the conduction band
mh = 0.57;					% effektive mass of the valence band
nue_e = 6;              % number of conduction band minima
mye_i = 1340;           % electron mobility in intrinsic material at 300 K
myh_i = 461;            % hole mobility in intrinsic material at 300 K

ED = EH * me/epsrel^2;  % donor binding energy in eV
aBe = aB0 * epsrel/me;  % donor electron radius in m
EA = EH * mh/epsrel^2;  % acceptor binding energy in eV
aBh = aB0 * epsrel/mh;  % acceptor hole radius in m

% calculation for 300 K
Nc_0 = 2 * nue_e * me^1.5 * 1.254708e19;  % effective density of conduction band states at 300 K
Nv_0 = 2 * mh^1.5 * 1.254708e19;          % effective density of valence band states at 300 K
ni_square_0 = Nc_0 * Nv_0 * exp(-EG/kT); 
ni_0 = sqrt(ni_square_0);                 % intrinsic carrier concentration at 300 K
