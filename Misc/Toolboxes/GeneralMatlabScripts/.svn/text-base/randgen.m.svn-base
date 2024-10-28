function Rand=randgen(Dist,N,ProbType,Method,NewSeed);
%---------------------------------------------------------------------------
% randgen function           General random numbers generator.
%                         Generate vector of random numbers from a given
%                         numerical distribution.
% Input  : - The distribution. This is a two column
%            vector in which the first column is the
%            independent variable and the second column is the
%            Probability Density Function (PDF) or CDF of the distribution.
%          - Number of elements in the random number vector.
%          - probability type:
%            'd' : differential PDF (default).
%            'c' : cumulative CDF.
%          - Interpolation type in the distribution matrix:
%            'linear'  - linear interpolation (default).
%            'nearest' - nearest neighbor interpolation.
%            'spline'  - cubic spline interpolation.
%            'cubic'   - cubic interpolation.
%          - NewSeed parameter: {'y' | 'n'}
%            'y' : take new seed each in each session.
%            'n' : do not take new seed (default).
% Output : - Vector of random numbers.
% Tested : Matlab 5.3
%     By : Eran O. Ofek        June 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
if (nargin==2),
   ProbType = 'd';
   Method   = 'linear';
   NewSeed  = 'n';
elseif (nargin==3),
   Method   = 'linear';
   NewSeed  = 'n';
elseif (nargin==4),
   NewSeed  = 'n';
elseif (nargin==5),
   % no default
else
   error('Illegal number of input arguments');
end

%MaxD = max(Dist(:,1));
%MinD = min(Dist(:,1));

%--- convert differential prob. to cumulative ---
switch ProbType
 case 'd'
    %Dist(:,2) = cumsum(Dist(:,2));
    Dist(:,2) = cumtrapz(Dist(:,1),Dist(:,2));
 case 'c'
    % do nothing
 otherwise
    error('Unknown probabiliyu type');
end


%--- normalize peak probability to 1 ---
Dist(:,2) = Dist(:,2)./max(Dist(:,2));

%--- check for monotoniticity ---
Diff = diff(Dist(:,2));
I = find(Diff==0);
%--- add epsilon to non monotonic values ---
Nel = length(Dist(:,1));
Eps = 1e-6./(N.*Nel);
I = I + 1;
Dist(I,2) = Eps.*[1:1:length(I)].' + Dist(I,2);

%--- normalize peak probability to 1 ---
Dist(:,2) = Dist(:,2)./max(Dist(:,2));



switch NewSeed
 case 'y'
    rand('state',sum(100*clock));
 case 'n'
    % do nothing
 otherwise
    error('Unknown NewSeed parameter {y | n}');
end

Y = rand(N,1);

Rand = interp1(Dist(:,2),Dist(:,1),Y);

