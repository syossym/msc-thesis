function Mid=binsear_f(Fun,Y,Range,Tol,varargin);
%---------------------------------------------------------------------------
% binsear_f function      Binary search on monotonic function
%                       (Given monotonic function y=f(x) and y, search for x).
% Input  : - Function [i.e., Y=Fun(X)].
%          - Y value to search
%          - X range to start the search in.
%          - Relative tolerance, default is 1e-3.
%          - Additional optional parameters of Fun
% Output : - X value corresponds to the input Y value.
% Tested : Matlab 6.5
%     By : Eran O. Ofek        Feb 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
if (nargin==3),
   Tol = 1e-3;
else
   % do nothing
end

Mid = mean(Range);
Y1  = feval(Fun,Range(1),varargin{:});
Y2  = feval(Fun,Range(2),varargin{:});
Ym  = feval(Fun,Mid,varargin{:});

if (Y2>Y1),
   % ascending function
   Type = 'a';
else
   % descending function
   Type = 'd';
end


while (diff(Range)>(Mid.*Tol)),
   switch Type
    case 'a'
       if (Y>Ym),
          Range = [Mid, Range(2)];
       else
          Range = [Range(1), Mid];
       end
    case 'd'
       if (Y>Ym),
          Range = [Range(1), Mid];
       else
          Range = [Mid, Range(2)];
       end
    otherwise
       error('Unknwon Type Option');
   end
   Mid = mean(Range);
   Y1  = feval(Fun,Range(1),varargin{:});
   Y2  = feval(Fun,Range(2),varargin{:});
   Ym  = feval(Fun,Mid,varargin{:});
end
