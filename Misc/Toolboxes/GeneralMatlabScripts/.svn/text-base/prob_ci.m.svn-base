function [MPX,SigX,OtherX]=prob_ci(ProbFun,AddProb,Eps);
%----------------------------------------------------------------------
% prob_ci function    Calculate 1,2,3\sigma (and others) confidence
%                   interval from a given numerical probability
%                   function.
% Input  : - Probability function [X,P]
%          - Additional probabilities values for which to calculate X.
%          - Epsilon value (for non monotonic series) default is 1e-12.
% Output : - Most probable X.
%          - [left sided 3\sigma, 2\sigma, 1\sigma, median,
%             right sided 1\sigma, 2\sigma, 3\sigma].
%          - X for the additional probability values.
% Tested : Matlab 5.3
%     By : Eran O. Ofek              October 2001
%                           last update: Feb 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------
if (nargin==1 | nargin==2),
   Eps = 1e-12;
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end

DefProb = [normcdf(-3,0,1), normcdf(-2,0,1), normcdf(-1,0,1),...
           normcdf( 0,0,1),... 
 	   normcdf( 1,0,1), normcdf( 2,0,1), normcdf( 3,0,1)];
% 0.0013    0.0228    0.1587    0.5000    0.8413    0.9772    0.9987

InterpMethod = 'linear';

% normalize probability function
X    = ProbFun(:,1);
P    = ProbFun(:,2)./trapz(X,ProbFun(:,2));
P    = (P + Eps);
P    = P./trapz(X,P);
CumP = cumtrapz(X,P);


[Val, Ind] = max(P);
MPX        = X(Ind);

SigX = interp1(CumP,X,DefProb,InterpMethod);

if (nargin>1),
   OtherX = interp1(CumP,X,AddProb,InterpMethod);
end
