function myh = myh(N,T)
% function myh(N,T)
% Löcherbeweglichkeit im Silizium
% Input:    N - Dichte der ionisierten Störstellen (ND+ + NA-)
%           T - Temperatur in Kelvin (Default: 300 K)       
% Output:   Beweglichkeit /u(T) in cm^2/Vs

% nach Shur M (1990) Physics of Semiconductor Devices. Prentice Hall, Englewood Cliffs, p. 76 
% und Roulston D J (1990) Bipolar Semiconductor Devices. Mc Graw Hill, New York, p. 7

if nargin < 2, T=300; end
N = N(:)'; T = T(:)';
Trel = T./300;
mym_h = 54.3 * (Trel).^-.57;
myo_h = 1.36e8 * T.^-2.23;
Ncp = 2.35e17 * Trel.^2.4;
ny = .88 * Trel.^-.146; 
myh = mym_h + myo_h./(1 + (N./Ncp).*ny);