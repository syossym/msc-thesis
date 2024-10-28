function X=mode(Vec,Method,Nbin,FitRange);
%--------------------------------------------------------------------------
% mode function       Estimate the mode value in a given vector
%                   (The most probable value).
% Input  : - Vector.
%          - Method :
%            'hist' - calculate mode from peak of histogram.
%            'fit'  - calculate mode from parabola fit to histogram.
%          - number of bins used in mode estimate,
%            (default is: floor(number of elements/10)+1).
%          - Fit range [Xlow, Xhigh]
% Output : - The mode value.
% Tested : Matlab 5.3
%     By : Eran O. Ofek               December 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
N = length(Vec);

if (nargin==2),
   Nbin = floor(N./10)+1;
   FitRange = [NaN, NaN];
elseif (nargin==3),
   FitRange = [NaN, NaN];
elseif (nargin==4),
   % no default
else
   error('Illigal number of input arguments');
end

switch Method
 case 'hist'
    Min = min(Vec);
    Max = max(Vec);
    [Xhist,Nhist]=realhist(Vec,[Min, Max, Nbin],'n');
    [Val, Ind] = max(Nhist);
    X = Xhist(Ind);
 case 'fit'
    Min = min(Vec);
    Max = max(Vec);
    [Xhist,Nhist] = realhist(Vec,[Min, Max, Nbin]);
    I             = find(Xhist>=FitRange(1) & Xhist<=FitRange(2));
    [Par,ParErr]  = fitpoly(Xhist(I),Nhist(I),sqrt(Nhist(I)),2);
    X             = -Par(2)./(2.*Par(3));
 otherwise
    error('Unknown Method option');
end
 
