function [Range]=rangend(Data);
%----------------------------------------------------------------------
% rangend function         Return the global Range of a N-D matrix
% Input  : - N-D matrix
% Output : - Value of global Range.
% Tested : Matlab 7.0
%     By : Eran O. Ofek      July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
%----------------------------------------------------------------------

SizeData = size(Data);
Vec = reshape(Data, [prod(SizeData), 1]);
[Range] = range(Vec);

