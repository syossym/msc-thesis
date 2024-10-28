function [MatchList,UnMatch1,UnMatch2]=match_lists_time(List1,List2,Thresh,Columns)
%--------------------------------------------------------------------------
% match_lists_time function     Match lists according to 1-D position
%                                 (e.g., time).
% Input  : - First list to match
%          - Second list to match
%          - position threshold [in 1-D position units]
%          - Columns to match, default is [1, 1] position matching with
%            the first column (in both lists).
% Output : - Matched list, in each raw, the columns of List1 are folowed
%            by the columns of List2, and than by the
%            1-D position distance between points, and the serial index,
%            to distinguish between points with more than one counterpart,
%            the point in List2 that is closest to the point in List1 get
%            the index 1 and so on.
%          - Unmatched raws in List1
%          - Unmatched raws in List2
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      January 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
if (nargin==3),
   Columns = [1,1];
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (length(Columns)==1),
   Columns = [Columns, Columns];
end

ColT1 = Columns(1);
ColT2 = Columns(2);

N1 = size(List1,1);
N2 = size(List2,1);


Jun1 = 0;
J    = 0;
MatchedJun2 = zeros(N2,1);
for I=1:1:N1,

   Dist = List1(I,ColT1)-List2(:,ColT2);

   Imatch = find(abs(Dist)<=Thresh);
   Nmatch = length(Imatch);
   %if (Nmatch>1),
   %   error('More than one raw in List2 is matched to List1');
   %end

   if (isempty(Imatch)==1),
      % no match
      Jun1 = Jun1 + 1;
      UnMatch1(Jun1,:) = List1(I,:);
   else
      % match
      J = J + 1;
      % sort Imatch by distance order
      [SortedDist, SortedDistInd] = sort(Dist(Imatch));
      MatchList(J:J+Nmatch-1,:) = [ones(Nmatch,1)*List1(I,:), List2(Imatch(SortedDistInd),:),SortedDist,[1:1:Nmatch].'];
      MatchedJun2(Imatch) = 1;
      J = J + Nmatch - 1;
   end
end

if (nargout>2),
   Jun2 = find(MatchedJun2==0);
   UnMatch2 = List2(Jun2,:);
end
