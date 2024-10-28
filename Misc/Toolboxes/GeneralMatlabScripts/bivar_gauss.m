function F=bivar_gauss(X,Y,Par);
%------------------------------------------------------------------
% bivar_gauss function     Return the value of a normalized
%                        bivariate gaussian in a list of points.
% Input  : - X
%          - Y
%          - Parameters: [X0, Y0, SigmaX, SigmaY, Rho]
% Output : - Value of bivariate gaussian at [X,Y].
% Tested : Matlab 7.0
%     By : Eran O. Ofek              May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------

X0     = Par(:,1);
Y0     = Par(:,2);
SigmaX = Par(:,3);
SigmaY = Par(:,4);
Rho    = Par(:,5);

F = 1./(2.*pi.*SigmaX.*SigmaY.*sqrt(1-Rho.^2)) .* ...
    exp(-1./(2.*(1-Rho.^2)) .* ...
 	                       ((X-X0).^2./SigmaX.^2 + ...
                                (Y-Y0).^2./SigmaY.^2 - ...
				2.*Rho.*(X-X0).*(Y-Y0)./(SigmaX.*SigmaY)));


