function r=rotm(p,n);
%--------------------------------------------------------------------
% rotm function     rotation matrix about the X/Y/Z axis. 
% input  : - parameter. angle to rotate in radians.
%          - rotation matrix number 1, 2 or 3 only,
%            1 - for rotation about the X axis.
%            2 - for rotation about the Y axis.
%            3 - for rotation about the Z axis.
% output : - rotation matrix number 1.
%    By  Eran O. Ofek           February 1994
%--------------------------------------------------------------------
if nargin~=2,
   error('2 arg only');
end
% P(p)
if n==1,
   r = [1, 0, 0; 0, cos(p), -sin(p); 0, sin(p), cos(p)];
% Q(p)
elseif n==2,
   r = [cos(p), 0, sin(p); 0, 1, 0; -sin(p), 0, cos(p)];
% R(p)
elseif n==3,
   r = [cos(p), -sin(p), 0; sin(p), cos(p), 0; 0, 0, 1];
else
   error('n must be 1, 2 or 3 only');
end

   
