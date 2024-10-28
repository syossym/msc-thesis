function DV=differr(V1,V2);
%--------------------------------------------------------------------
% differr function       Calculate differences and errors.
%                      If two input arguments V1, V2 are given then
%                      the program returns:
%                      [V1(:,1)-V2(:,1), sqrt(V1(:,2).^2+V2(:,2).^2)]
%                      If one input argument is given the the program
%                      differentiate successive values.
% Input  : - Two column data matrix [Val, Err].
%          - Two column data matrix [Val, Err].
%            (optional).  
% Output : - Two column matrix in which the first column
%            is the difference between the 'values' in 
%            the two matrices
%            and the second column is the error as calculated
%            from the two 'values errors'.
% Tested : Matlab 5.3
%     By : Eran O. Ofek           Febuary 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
ValCol = 1;
ErrCol = 2;
if (nargin==1),
   N = length(V1(:,ValCol));
   DV = zeros(N-1,2);
   DV(:,1) = V1(1:N-1,ValCol)-V1(2:N,ValCol);
   DV(:,2) = sqrt(V1(1:N-1,ErrCol).^2 + V1(2:N,ErrCol).^2);
elseif (nargin==2),
   N = length(V1(:,ValCol));
   DV = zeros(N,2);
   DV(:,1) = V1(:,ValCol) - V2(:,ValCol);
   DV(:,2) = sqrt(V1(:,ErrCol).^2 + V2(:,ErrCol).^2);
else
   error('Illigal number of input arguments');
end
