function [Max,MaxIndVec]=maxnd(Data);
%----------------------------------------------------------------------
% maxnd function          Return the global maximum of a N-D matrix
%                       and its indices.
% Input  : - N-D matrix
% Output : - Value of global Maximum
%          - Vector that contains the indices of the global maximum
%            in the N-D matrix.
% Tested : Matlab 7.0
%     By : Eran O. Ofek      July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
%----------------------------------------------------------------------

SizeData = size(Data);
Vec = reshape(Data, [prod(SizeData), 1]);
[Max, MaxInd] = max(Vec);

ND = length(SizeData);

StrCom = '[I1 ';
for Idim=2:1:ND,
   StrTemp = sprintf('I%d ',Idim);
   StrCom  = [StrCom,',',StrTemp];
end
StrInd = [StrCom,']'];
StrCom = [StrInd,'=ind2sub(SizeData,MaxInd);'];
eval(StrCom);
MaxIndVec=eval(StrInd);

%[I1, I2, I3] = ind2sub(SizeData,MaxInd)


   
