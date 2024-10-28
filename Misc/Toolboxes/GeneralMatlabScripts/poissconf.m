function CI=poissconf(N,CL);
%---------------------------------------------------------------------------
% poissconf function     Given the number of observed events, assuming
%                      Poisson statistics, calculate the two sided upper
%                      and lower confidence intervals.
% Input  : - Coloumn vector of numbers.
%          - Two sided confidence limit, default is the 1\sigma 0.6827 CL.
%            Can get also '1' - for 1\sigma 0.6827 CL.
%                         '2' - for 2\sigma 0.9545 CL.
%                         '3' - for 3\sigma 0.9973 CL.
% Output : - Two column matrix of lower and upper confidence intervals.
%            Left column for lower CI, right column for right CI. 
% Examples : CI=poissconf([0;1;2;3;4]);
% Reference: Gehrels, N. 1986, ApJ 303, 336
% Tested : Matlab 5.3
%     By : Eran O. Ofek         December 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
%error('bugs');
SIGMA1 = 0.6827;
SIGMA2 = 0.9545;
SIGMA3 = 0.9973;

LSTEP     = 1e-8;
THRESHOLD = 1e-8;

PoissFun = inline('sum((L.^X) .* exp(-L)./gamma(X+1))','L','X'); 

LenN = length(N);

if (nargin==1),
   CL = SIGMA1;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isstr(CL)==1),
   switch CL
    case '1'
       CL = SIGMA1;
    case '2'
       CL = SIGMA2;
    case '3'
       CL = SIGMA3;
    otherwise
       error('poissconf: Unkown CL option');
   end
end

if (CL<0 | CL>1),
   error('CL should be in the range 0 to 1');
end


% convert two sided CL to singl sided CL
CL = 1 - 0.5.*(1-CL);

CI = zeros(LenN,2);
for I=1:1:LenN,
   %--------------
   % lower CI (LL)
   %--------------
   if (N(I)==0),
      CI(I,1) = NaN;
   else
      X = [0:1:N(I)-1].';
   
      LL0 = N(I)-0.25;
      LL1 = LL0 + LSTEP;

      % Solve using Newton-Rapson method
      DelL = 1 + THRESHOLD;
      while (abs(DelL)>THRESHOLD),
         F0 = PoissFun(LL0,X);
         F1 = PoissFun(LL1,X);

         Df = (F1 - F0)./LSTEP;

         DelL = -(F0 - CL)./Df;

         LL0 = LL0 + DelL;
         LL1 = LL0 + LSTEP;
      end
      CI(I,1) = LL0;
   end

   %--------------
   % upper CI (LU)
   %--------------
   X = [0:1:N(I)].';
   
   LU0 = N(I);
   LU1 = LU0 + LSTEP;

   % Solve using Newton-Rapson method
   DelL = 1 + THRESHOLD;
   while (abs(DelL)>THRESHOLD),
      F0 = PoissFun(LU0,X);
      F1 = PoissFun(LU1,X);

      Df = -(F1 - F0)./LSTEP;

      DelL = (F0 - (1-CL))./Df;

      LU0 = LU0 + DelL;
      LU1 = LU0 + LSTEP;
   end
   CI(I,2) = LU0;

end

