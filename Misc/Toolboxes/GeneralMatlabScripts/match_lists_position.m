function [MatchList,UnMatch1,UnMatch2]=match_lists_position(List1,List2,Thresh,Columns,CosD)
%--------------------------------------------------------------------------
% match_lists_position function     Match lists according to 2-D position.
% Input  : - First list to match
%          - Second list to match
%          - position threshold [in position units]
%          - Columns to match, default is position matching with
%            the columns [1 2] (in both lists).
%            In case four numbers are given, its taken as the columns in List1
%            and in List2, respectively.
%          - use spherical distance : {'y' | 'n'} (default is 'y').
%            Assuming position is in radians, and that first column
%            is longitude and second column is latitude.
% Output : - Matched list, in each raw, the columns of List1 are folowed
%            by the columns of List2, and than by the
%            Euclidean distance between points, and the serial index,
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
   Columns = [1,2,1,2];
   CosD = 'y';
elseif (nargin==4),
   CosD = 'y';
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end




if (length(Columns)==2),
   Columns = [Columns, Columns];
end

ColX1 = Columns(1);
ColY1 = Columns(2);
ColX2 = Columns(3);
ColY2 = Columns(4);

N1 = size(List1,1);
N2 = size(List2,1);


Jun1 = 0;
J    = 0;
MatchedJun2 = zeros(N2,1);
for I=1:1:N1,
   if (CosD=='y'),
      Dist = sphere_dist(List1(I,ColX1), List1(I,ColY1), List2(I,ColX2), List2(I,ColY2));
   else
      Dist = sqrt((List1(I,ColX1)-List2(:,ColX2)).^2 + (List1(I,ColY1)-List2(:,ColY2)).^2);
   end

   Imatch = find(Dist<=Thresh);

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
