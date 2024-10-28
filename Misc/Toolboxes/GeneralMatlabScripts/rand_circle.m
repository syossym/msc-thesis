function [X,Y,R]=rand_circle(N);
%---------------------------------------------------------------------------
% rand_circle function         Generate random number equally distributed
%                            inside a unit circle.
% Input  : - Number of random-numbers to generate.
% Output : - X.
%          - Y.
%          - R (radius).
% Tested : MATLAB 5.1
%     By : Eran O. Ofek         Feb 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------

Ngen = ceil(N.*1.3 + sqrt(5.*N) + 20);

rand('state',sum(100*clock));

X = rand(Ngen,1).*2 - 1;
Y = rand(Ngen,1).*2 - 1;

R = sqrt(X.^2 + Y.^2);
I = find(R<=1);

X = X(I(1:N));
Y = Y(I(1:N));
R = R(I(1:N));
