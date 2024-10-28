function Prob=prob2find_inr(Density,R);
%----------------------------------------------------------------------------
% prob2find_inr function     Given a density (number per unit area),
%                          and a distance from a point, calculate the
%                          the probability to find an object within radius R
%                          from the point (assuming Poisson statistics).
% Input  : - Object density [Nuber/Area]
%          - Radius [Dist]
% Output : - Probability to find an object within radius R.
% Tested : Matlab 5.3
%     By : Eran O. Ofek           Feb 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------
Prob = 1-exp(-pi.*Density.*R.^2);
