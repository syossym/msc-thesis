function [M,E]=wmean(Vec)
%--------------------------------------------------------------------
% wmean function     Wighted Mean     (Ignoring NaNs).
% input  : - Vector in which the first column is the data, and
%            the second column is the std. errors.
% output : - Wighted Mean
%	   - Error.
%    By  Eran O. Ofek           June 1998
%--------------------------------------------------------------------


C_D = 1;
C_E = 2;

I = find(isnan(Vec(:,C_D))==0 & isnan(Vec(:,C_E))==0);

E = sqrt(1./sum(1./Vec(I,C_E).^2));
M = sum(Vec(I,C_D)./(Vec(I,C_E).^2))./sum(1./(Vec(I,C_E).^2));
