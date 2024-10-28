function SN=sn_calc(ExpTime,StarMag,varargin);
%---------------------------------------------------------------------------
% sn_calc function       Calculate eproximate Signal to Noise (S/N) ratio.
%                      Given telescope parameters, calculate the S/N of
%                      a point source.
% Input  : - Exposure time [sec].
%          - Magnitude.
%          - Pairs of keywords and value for each keyword.
%            Keyword, one of the following keywords:
%            'TA' - Telescope aperture (assuming unfiltered
%                   telescope with CCD QE~80% - Wise case).
%                   Default 1.0 m.
%            'ZP' - Telescope/imager zero point (override 'TA').
%                   Default 22.2 [mag adu/sec]
%            'SB' - Sky brightness.
%                   Default 20.0 [mag adu/sec]
%            'Ga' - CCD Gain.
%                   Default 8.42 [electron/adu]
%            'RN' - CCD read out noise.
%                   Default 6.5 [electron]
%            'Sc' - CCD scale [arcsec/pixel].
%                   Default 0.7 [arcsec/pixel]
%            'Se' - Seeing [arcsec].
%                   Default 2.3
%            'RP' - photometry aperture radius [pixel].
%                   One seeing FWHM [pixel].
%            'AM' - Air mass.
%                   Default 1.0
%            'EC' - Extinction coef.
%                   Default 0.3
%            'De' - Dark current [electron/sec/pixel].
%                   Default 0
%            'DC' - Dark current [adu/sec/pixel] (override 'De').
%                   Default 0
%            'ND' - Number of dark frames.
%                   Default 1
%            'BA' - Background area [pixel].
%                   Default 300
%          - The value for previous keyword.
% Output : - S/N
% Example: SN=sn_calc(100,[19;20],'Se',1.0,'AM',2.0);
% Tested : MATLAB 5.3
%     By : Eran O. Ofek                September 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------

% default values
TA = 1.0;
ZP = 22.2;
SB = 20.0;
Ga = 8.42;
RN = 6.5;
Sc = 0.7;
Se = 2.3;
RP = Se./Sc;
AM = 1.0;
EC = 0.3;
De = 0.0;
DC = 0.0;
ND = 1;
BA = 300;


Narg = length(varargin);
if (floor(0.5.*Narg)==(0.5.*Narg)),
   % OK
else
   error('Should give one value per keyword');
end

ZPo = 0;
DCo = 0;
RPo = 0;
for I=1:2:Narg-1,
   switch varargin{I}
    case 'TA'
       TA = varargin{I+1};
    case 'ZP'
       ZP = varargin{I+1};
       ZPo= 1;
    case 'SB'
       SB = varargin{I+1};
    case 'Ga'
       Ga = varargin{I+1};
    case 'RN'
       RN = varargin{I+1};
    case 'Sc'
       Sc = varargin{I+1};
    case 'Se'
       Se = varargin{I+1};
    case 'RP'
       RP = varargin{I+1};
       RPo= 1;
    case 'AM'
       AM = varargin{I+1};
    case 'EC'
       EC = varargin{I+1};
    case 'De'
       De = varargin{I+1};
    case 'DC'
       DC = varargin{I+1};
       DCo= 1;
    case 'ND'
       ND = varargin{I+1};
    case 'BA'
       BA = varargin{I+1};
    otherwise
       error('Unknown keyword');
   end
end


SkyB        = SB;
if (ZPo==1),
   ZeroP    = ZP;
else
   ZeroP    = 22.20 + 2.5.*log10(TA.*TA);   % mag->adu/sec
end
Gain         = Ga;
ReadNoise    = RN;   % electrons
if (DCo==1),
   AduSec_Dark  = DC;   % adu/sec/pixel
   ElecSec_Dark = AduSec_Dark.*Gain;   % electron/sec/pixel
else
   AduSec_Dark  = De./Gain;   % adu/sec/pixel
   ElecSec_Dark = De;   % electron/sec/pixel
end
Scale        = Sc; %[0.5:0.05:6]';   %0.696;  % arcsec/pixel
Seeing       = Se;  % arcsec
FWHM         = Seeing./Scale;   % FWHM in units of pixels
if (RPo==1),
   R_optimal = RP;
else
   R_optimal = FWHM;
end   
AirMass      = AM;
ExtinCoeff   = EC;
N_Dark_Frame = ND;
BackAreaPix  = BA;


%convert ZeroP to electrons:
ZeroP = ZeroP + 2.5.*log10(Gain);

Sigma = FWHM./(2.*1.1774);  % in pixels units

AperAreaPix  = pi.*R_optimal.^2;
Is = find(AperAreaPix<1);
Il = find(AperAreaPix>=1);
AperAreaPix(Is) = 1;
R_optimal(Is)   = sqrt(1./pi);  %0.5642

%StarLightFracInAper = (1 - exp(-R_optimal.^2./(2.*Sigma.^2)));
StarLightFracInAper = (1 - exp(-4.*log(2).*(R_optimal./FWHM).^2));



RelStarMag = ZeroP - (StarMag + ExtinCoeff.*(AirMass-1));
RelSkyB    = ZeroP - SkyB;

ElecSec_Sky  = 10.^(0.4.*RelSkyB);  % elec/sec/sq. arcsec.
ElecSec_Sky  = ElecSec_Sky.*Scale.*Scale;
ElecSec_Star = 10.^(0.4.*RelStarMag);
ElecSec_Star = ElecSec_Star.*StarLightFracInAper;

% signal
Signal = ElecSec_Star.*ExpTime;

% noise
Total_Dark = (ElecSec_Dark + ElecSec_Dark./sqrt(N_Dark_Frame)).*ExpTime.*AperAreaPix;
Total_Sky  = ElecSec_Sky.*ExpTime.*AperAreaPix;
Noise  = sqrt(Signal + (1+AperAreaPix./BackAreaPix).*(Total_Sky + Total_Dark + AperAreaPix.*ReadNoise.^2));

SN = Signal./Noise;




