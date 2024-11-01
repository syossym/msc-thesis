function OutE=convert_energy(InE,InUnit,OutUnit);
%---------------------------------------------------------------
% convert_energy function    Convert between different energy
%                          units.
% Input  : - Energy.
%          - Input system:
%            'erg'   - ergs
%            'J'     - Jouls
%            'Hz'    - Frequency [1/s]
%            'A'     - Wavelength [Ang]
%            'eV'    - Electron volts [h nu/q]
%            'T'     - Temperature [K]
%            'me'    - Electron mass [E/m_e]
%            'mp'    - Proton mass [E/m_p]
%            'cal'   - calorie (4.184 J)
%            'Btu'   - (1.055x10^3 J)
%            'kWh'   - kilowatt-hour (3.6x10^6 J)
%            'TNT'   - one ton of TNT (4.2x10^9 J) 
%            'gr'    - Energy equivalent of 1 gram of matter (9x10^13 J)
%          - Output system (see input system for options).
% Output : - Energy in output system.
% Tested : Matlab 7.0
%     By : Eran O. Ofek            Feb 2006
%                                  Aug 2006 (bug correction)
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------
Erg2Erg = 1;
Erg2J   = 1e-7;
Erg2Hz  = 1.5092e26;
Erg2A   = 1.9864e-8; 
Erg2eV  = 6.2415e11;
Erg2T   = 7.2430e15;
Erg2me  = 1.2214e6;
Erg2mp  = 665.214577;
Erg2cal = Erg2J./4.184;
Erg2Btu = Erg2J./1.055e3;
Erg2kWh = Erg2J./3.6e6;
Erg2TNT = Erg2J./4.2e9;
Erg2gr  = get_constant('c','cgs').^-2;

Relation = 'lin';

switch InUnit
 case 'erg'
    ConvFactor = Erg2Erg;
 case 'J'
    ConvFactor = Erg2J;
 case 'Hz'
    ConvFactor = Erg2Hz;
 case 'A'
    Relation   = 'inv';
    ConvFactor = Erg2A;
 case 'eV'
    ConvFactor = Erg2eV;
 case 'T'
    ConvFactor = Erg2T;
 case 'me'
    ConvFactor = Erg2me;
 case 'mp'
    ConvFactor = Erg2mp;
 case 'cal'
    ConvFactor = Erg2cal;
 case 'Btu'
    ConvFactor = Erg2Btu;
 case 'kWh'
    ConvFactor = Erg2kWh;
 case 'TNT'
    ConvFactor = Erg2TNT;
 case 'gr'
    ConvFactor = Erg2gr;
 otherwise
    error('Unknown InUnit Option');
end

switch Relation
 case 'lin'
    ErgE = InE./ConvFactor;
 case 'inv'
    ErgE = ConvFactor./InE;
 otherwise
    error('Unknown Relation Option');
end


Relation = 'lin';
switch OutUnit
 case 'erg'
    ConvFactor = Erg2Erg;
 case 'J'
    ConvFactor = Erg2J;
 case 'Hz'
    ConvFactor = Erg2Hz;
 case 'A'
    Relation   = 'inv';
    ConvFactor = Erg2A;
 case 'eV'
    ConvFactor = Erg2eV;
 case 'T'
    ConvFactor = Erg2T;
 case 'me'
    ConvFactor = Erg2me;
 case 'mp'
    ConvFactor = Erg2mp;
 case 'cal'
    ConvFactor = Erg2cal;
 case 'Btu'
    ConvFactor = Erg2Btu;
 case 'kWh'
    ConvFactor = Erg2kWh;
 case 'TNT'
    ConvFactor = Erg2TNT;
 case 'gr'
    ConvFactor = Erg2gr;
 otherwise
    error('Unknown InUnit Option');
end


switch Relation
 case 'lin'
    OutE = ErgE.*ConvFactor;
 case 'inv'
    OutE = ConvFactor./ErgE;
 otherwise
    error('Unknown Relation Option');
end
