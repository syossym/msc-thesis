function FS=fresnels(T);
%--------------------------------------------------------
% fresnels function      return: sin(0.5*pi*T^2)
% Input  : - Parameter.
% Output : - Fresnel sine function.
% Tested : Matlab 5.3
%     By : Eran O. Ofek       April 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------
FS = sin(0.5.*pi.*T.^2);
