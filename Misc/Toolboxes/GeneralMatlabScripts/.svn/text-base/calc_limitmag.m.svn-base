function LimitMag=calc_limitmag(ExpTime,LimitSN,varargin);
%---------------------------------------------------------------------------
% sn_calc function       Calculate eproximate limiting magnitude.
%                      Given telescope parameters, and S/N calculate the
%                      limiting magnitude for a point source.
% Input  : - Exposure time [sec].
%          - Signal to Noise.
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
% See also: calc_sn.m
% Tested : MATLAB 5.3
%     By : Eran O. Ofek                May 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
Dmag      = 0.01;
StartMag  = 20;
SN_Thresh = 0.01;


%--- solve using Newton-Rapson method ---
Mag0 = StartMag;
SN0  = LimitSN+10.*SN_Thresh;

while (abs(LimitSN-SN0)>SN_Thresh),
   SN0 = sn_calc(ExpTime,Mag0,varargin{:});

   Mag1 = Mag0 + Dmag;
SN1 = sn_calc(ExpTime,Mag1,varargin{:});

   DSN_Dmag = (Mag0 - Mag1)./(SN0 - SN1);

   Mag0 = Mag0 - DSN_Dmag.*(SN0 - LimitSN);


end

LimitMag = Mag0;
