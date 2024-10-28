function i=min_sear(x,v);
%--------------------------------------------------------------------
% min_sear function        searching for value position in a given
%                        vector by minimazation. 
% input  : - vector.
%          - value to search.
% output : - index of closest value.
%    By  Eran O. Ofek           September 1994
%--------------------------------------------------------------------
N = length(x);
if nargin~=2,
   error('2 args only');
end
[m,i] = min(abs(x-v));
