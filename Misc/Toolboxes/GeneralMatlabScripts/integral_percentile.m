function [LeftLimit,RightLimit,TotInt]=integral_percentile(X,Y,Per,Type,Method);
%-------------------------------------------------------------------------
% integral_percentile function   Given a function (X,Y), calculate
%                              the limits of the integal which contains
%                              a given percentile of the total integral.
% Input  : - X
%          - Y
%          - Vector of percentile, for each to calculate the limits
%            of the integral which contains that percentiles of the total
%            integral.
%          - Percentile type:
%            'c'  - central (default). The limits are such that the integarl
%                   from min(X) to the left limit is equal to the
%                   integral from the right limit to max(X).
%            'l'  - left. The limits are such that the left limit is 0.
%            'r'  - left. The limits are such that the right limit is 0.
%          - Method:
%            'Interp'  - work well for smooth functions.
%            'First'   - 
% Output : - Vector of left limits.
%          - Vector of right limits.
%          - Total integral.
% Tested : Matlab 7.0
%     By : Eran O. Ofek          December 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------
InterpMethod = 'linear';

DefType   = 'c';
DefMethod = 'First';
if (nargin==3),
   Type   = DefType;
   Method = DefMethod;
elseif (nargin==4),
   Method = DefMethod;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end



CumInt    = cumtrapz(X,Y);
TotInt    = max(CumInt);
PerCumInt = CumInt./TotInt;

switch Type
 case 'c'
    switch Method
     case 'Interp'
        LeftLimit   = interp1(PerCumInt,X,    0.5.*(1-Per), InterpMethod);
        RightLimit  = interp1(PerCumInt,X,1 - 0.5.*(1-Per), InterpMethod);
     case 'First'
        for Ip=1:1:length(Per),
           I = find( (PerCumInt - 0.5.*(1-Per(Ip))) > 0);
           LeftLimit(Ip,1) = X(I(1));
           I = find( (PerCumInt - (1-0.5.*(1-Per(Ip)))) < 0);
           RightLimit(Ip,1) = X(I(end));
        end
     otherwise
        error('Unknown Method Option');
    end
 case 'l'
    RightLimit = interp1(PerCumInt,X,Per, InterpMethod);
    LeftLimit  = zeros(size(RightLimit));
 case 'r'
    LeftLimit  = interp1(PerCumInt,X,1-Per, InterpMethod);
    RightLimit = ones(size(LeftLimit)).*max(X);
 otherwise
    error('Unknown Type Option');
end
