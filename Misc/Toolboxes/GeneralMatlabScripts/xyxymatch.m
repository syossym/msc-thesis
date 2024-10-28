function [MatchInd,MatchInd2]=xyxymatch(List1,List2,ColCooList1,ColCooList2,DistThresh,Separation,CooType);
%----------------------------------------------------------------------
% xyxymatch function     Match two lists by 2-d coordinates (spherical
%                      or planner). The program indicate about possible
%                      confusion in the matching (i.e., one to many...).
% Input  : - First list.
%          - Second list -- sorted by Y coordinate! --
%          - Coloumns of X and Y coordinates in first list,
%            if [], then use default. Default is [1 2].
%          - Coloumns of X and Y coordinates in second list,
%            if [], then use default. Default is [1 2].
%          - Matching radius threshold in coordinates units.
%          - Minimum separation - The number of objects in List2
%            found within minimum-separation radius from each object
%            in List1 are indicated.
%          - Coordinates type:
%            'sphere' - Spherical coordinates. Coordinates should be
%                       given in radians.
%            'plane'  - plan coordinates. Default.
% Output : - Matched indices [Ind1, Ind2, Dist, PA, NumberTh, NumberSep, MultMatch],
%            where "Ind2" is the index of line in List2 that was
%            identified (closest) to the line "Ind1" in List1.
%              "NumberTh" is the number of objects in List2
%            found within threshold radius from each object
%            in List1.
%              "NumberSep" is the number of objects in List2
%            found within minimum-separation radius from each object
%            in List1.
%              "MultMatch" is a {0 | 1} flag indicating if the object
%            in List2 (i.e., Ind2) was matched to more than one object
%            in List1 (i.e., Ind1).
%              Objects in List1 that doesn't have a counterpart within
%            threshold radius are marked as [Ind1, NaN, NaN, NaN, 0 NumberSep 0].
%          - All the objects in List2 that doesn't have a counterpart in List1
%            are added to the MatchInd2 list as:
%            [Ind1, Ind2, Dist, PA, NumberTh, NumberSep, NaN]
% Tested : Matlab 7.0
%     By : Eran O. Ofek                September 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------
SortedCoo = 'Dec';
Col_Ind1      = 1;
Col_Ind2      = 2;
Col_Dist      = 3;
Col_PA        = 4;
Col_Nth       = 5;
Col_Nsep      = 6;
Col_MultMatch = 7;
Col_Tot       = 7;


if (nargin==6),
   CooType = 'plane';
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(ColCooList1)==1),
   ColCooList1 = [1 2];
end
if (isempty(ColCooList2)==1),
   ColCooList2 = [1 2];
end


N1 = size(List1,1);
N2 = size(List2,1);

FlagList2 = zeros(N2,1);       % Flag objects in List2: 0-not found in list1; 1-found in List1
MatchInd = zeros(N1,Col_Tot);  %[Ind1, Ind2, Dist, PA, NumberTh, NumberSep, MultMatch]
for I1=1:1:N1,

   [LinesSep,DistL,PAL] = cat_search(List2, ColCooList2, List1(I1,ColCooList1), Separation, 'circle', SortedCoo, CooType);
   NumberSep = length(LinesSep);
   if (NumberSep>0),
      LinesTh = LinesSep(find(DistL<DistThresh));   % indices of objects within threshold-radius
      [MinDist,MinInd] = min(DistL);
      I2               = LinesSep(MinInd);
      PA               = PAL(MinInd);
      FlagList2(I2)    = 1;
   else
      LinesTh = [];
      MinDist = NaN;
      I2      = NaN;
      PA      = NaN;
   end
   NumberTh = length(LinesTh);

   MatchInd(I1,:) = [I1, I2, MinDist, PA, NumberTh, NumberSep, NaN];
end



%--- populate 7th column in MatchInd (MultMatch) ---
%[Ind1, Ind2, Dist, PA, NumberTh, NumberSep, MultMatch]
%              "MultMatch" is a {0 | 1} flag indicating if the object
%            in List2 (i.e., Ind2) was matched to more than one object
%            in List1 (i.e., Ind1).
%SortedInd2 = sort(MatchInd(:,Col_Ind));
MatchInd(:,Col_MultMatch)       = 0;
for I1=1:1:N1,
   CurInd2  = MatchInd(I1,Col_Ind2);
   Imatched = find(MatchInd(:,Col_Ind2)==CurInd2);
   if (length(Imatched)>1),
      MatchInd(Imatched,Col_MultMatch) = 1;
   %else
   %   MatchInd(I1,Col_MultMatch)       = 0;
   end
end


%--- Adding objects in List2 that are not in List1 ---
if (nargout>1),
   I2not_in_List1 = find(FlagList2==0);
   N2not_in_List1 = length(I2not_in_List1);
   MatchInd2      = zeros(N2not_in_List1,Col_Tot);
   
   % populate MatchInd2
   if (length(I2not_in_List1)==0),
      % do nothing
   else
      MatchInd2(:,Col_Ind2) = I2not_in_List1;
   end

   for K2=1:1:N2not_in_List1,
      I2 = I2not_in_List1(K2);
      [LinesSep,DistL,PAL] = cat_search(List1, ColCooList1, List2(I2,ColCooList2), Separation, 'circle', SortedCoo, CooType);
      NumberSep = length(LinesSep);
      if (NumberSep>0),
         LinesTh = LinesSep(find(DistL<DistThresh));   % indices of objects within threshold-radius
         [MinDist,MinInd] = min(DistL);
         I1               = LinesSep(MinInd);
         PA               = PAL(MinInd);
         %FlagList2(I2)    = 1;
      else
         LinesTh = [];
         I1      = NaN;
      end
      NumberTh = length(LinesTh);
   
      MatchInd2(K2,:) = [I1, I2, MinDist, PA, NumberTh, NumberSep, NaN];
   end

end
