function [x,y,q]=equipot(m1,m2,a,s,n,z);
%--------------------------------------------------------------------
% equipot function      calculate two body equipotanials map.
% Input  : - m1 : first star mass in solar mass.
%          - m2 : second star mass in solar mass.
%          - a : the separation between the two stars. in meters.
%          - s : scaling the result from min to max in units of
%            the staller separation.
%          - n : number of point in each axis. default is 15.
%          - z : the surface to work on. default is z=0 (orbital plane)  
% Output : - grid of x coordinate.
%          - grid of y coordinate.
%          - matrix of potential defined by the x/y grids.
% Tested : Matlab 5.3
%     By : Eran O. Ofek           May 1994
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
if nargin==4,
   z = 0;
   n = 15;
elseif nargin==5,
   z = 0;
elseif nargin>6,
   error('4, 5 or 6 args only');
elseif nargin<4,
   error('4, 5 or 6 args only');
end
G = 6.672e-11;
m1 = m1.*1.9891e30;
m2 = m2.*1.9891e30;
mu = m2./(m1 + m2);
q  = m2./m2;
om2= G.*(m1 + m2)./(a.*a.*a);
d  = s.*a./n;
x  = -s.*a:d:2.*s.*a;
y  = -s.*a:d:2.*s.*a;
for i=1:length(x),
   for j=1:length(y),
      t1 = G.*m1./sqrt(x(i).*x(i) + y(j).*y(j) + z.*z);
      t2 = G.*m2./sqrt((x(i) - a).*(x(i) - a) + y(j).*y(j) + z.*z);
      t3 = om2.*((x(i) - mu.*a).*(x(i) - mu.*a) + y(j).*y(j));
      q(i,j) = -t1 - t2 - 0.5.*t3;
   end
end
