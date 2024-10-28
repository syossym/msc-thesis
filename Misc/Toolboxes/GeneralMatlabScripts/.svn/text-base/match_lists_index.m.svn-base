function [MatchList,UnMatch1,UnMatch2]=match_lists_index(List1,List2,Columns)
%--------------------------------------------------------------------------
% match_lists_index function     Match lists according to index.
% Input  : - First list to match
%          - Second list to match
%          - Columns to match, default is index matching with
%            the first columns.
%            In case two numbers are given, its taken as the column in List1
%            and in List2, respectively.
% Output : - Matched list, in each raw, the columns of List1 are folowed
%            by the columns of List2.
%          - Unmatched raws in List1
%          - Unmatched raws in List2
% Example: [MatchList,UnMatch1,UnMatch2]=match_lists_index(List1,List2,[1 1]);
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      December 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
if (nargin==2),
   Columns = [1,1];
end

if (length(Columns)==1),
   Columns = [Columns, Columsn];
end

Col1 = Columns(1);
Col2 = Columns(2);

N1 = size(List1,1);
N2 = size(List2,1);


Jun1 = 0;
J    = 0;
MatchedJun2 = zeros(N2,1);
for I=1:1:N1,
   Imatch = find(List1(I,Col1)==List2(:,Col2));
   if (length(Imatch)>1),
      error('More than one raw in List2 is matched to List1');
   end

   if (isempty(Imatch)==1),
      % no match
      Jun1 = Jun1 + 1;
      UnMatch1(Jun1,:) = List1(I,:);
   else
      % match
      J = J + 1;
      MatchList(J,:) = [List1(I,:), List2(Imatch,:)];
      MatchedJun2(Imatch) = 1;
   end
end

if (nargout>2),
   Jun2 = find(MatchedJun2==0);
   UnMatch2 = List2(Jun2,:);
end
