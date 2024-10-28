function [NewData,Nboundary]=ensemble_from_error(Data,Errors,Boundary);
%--------------------------------------------------------------------------
% ensemble_from_error function    Generate normal distributed random
%                               ensamble realization given a mean value
%                               and the one (or two) sided standard
%                               deviation.
% Input  : - Data vector.
%          - One sigma error vector.
%            If two column matrix is given, than the
%            first column is for the lower-side errors, and the second
%            column is for the upper-side errors.
%          - Vector of boundary values [down, up], default is [-Inf, Inf].
%            if value is outside boundary than it is replaced by
%            the nearest boundary value.
% Output : - Vector of new data, normaly distributed around the original
%            data pints.
%          - Number of boundary replacments.
% Tested : Matlab 5.3
%     By : Eran O. Ofek
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
if (nargin==2),
   Boundary = [-Inf, Inf];
elseif (nargin==3),
   % no default
else
   error('Illegal number of input arguments');
end

if (length(Errors(1,:))==1),
   Errors = [Errors, Errors];
end

rand('state',sum(100.*clock));

N = length(Data);

Rand = randn(N,1);

IRandN = find(Rand<0);
IRandP = find(Rand>=0);

NewData = zeros(size(Data));

if (isempty(IRandN)==0),
   NewData(IRandN) = Data(IRandN) + Errors(IRandN,1).*Rand(IRandN);
end
if (isempty(IRandP)==0),
   NewData(IRandP) = Data(IRandP) + Errors(IRandP,2).*Rand(IRandP);
end

if (Boundary(1)==-Inf & Boundary(2)==Inf),
   % no boundary
   Nboundary = 0;
else
   Ib_l = find(NewData<Boundary(1));
   Ib_u = find(NewData>Boundary(2));
   if (isempty(Ib_l)==1),
      % no point is out of lower boundary
      Nl = 0;
   else
      NewData(Ib_l) = Boundary(1);
      Nl = length(Ib_l);
   end
   if (isempty(Ib_u)==1),
      % no point is out of upper boundary
      Nu = 0;
   else
      NewData(Ib_u) = Boundary(2);
      Nu = length(Ib_u);
   end
   Nboundary = Nl + Nu;
end
