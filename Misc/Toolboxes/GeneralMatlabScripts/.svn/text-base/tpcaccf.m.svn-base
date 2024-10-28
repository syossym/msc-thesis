function Cf=tpcaccf(List1,List2,TetaMax,DelTeta);
%--------------------------------------------------------------------
% tpccf function    two point comulative angular
%                 crosscorrelation function.
% input  : - first series matrix, first column contains
%            celecstial longitude, second column contains
%            celestial latitude and third column contain
%            error of position. All values in radians.
%            The first series should be the smaller list!
%          - second series matrix, first column contains
%            celecstial longitude, second column contains
%            celestial latitude and third column contain
%            error of position. All values in radians.
%          - TetaMax is the maximum angle to which calculate
%            the two point ccf. Given in radians.
%          - DelTeta is the interval in Teta. Given in radians.
% output : - Cros correlation funtion as function of Teta.
%            Two column matrix in which the first column is Teta
%            and the second for the CCF.
% reference :
%    By  Eran O. Ofek           April 1997
%--------------------------------------------------------------------

%CEr = 3;    % column number of position error
CRA = 1;    % column number of R.A.
CDec = 2;   % column number of Dec.

N1 = length(List1(:,CRA));    % number of elements in first list
N2 = length(List2(:,CRA));    % number of elements in second list


Cf = zeros(floor(TetaMax./DelTeta),2);    % initialize the correlation coef.

j = 0;
for Teta=DelTeta:DelTeta:TetaMax,
   j = j + 1;
   for i=1:1:N1,
      % find distance between objects
      Dist = acos(sin(List1(i,CDec)).*sin(List2(:,CDec))+cos(List1(i,CDec)).*cos(List2(:,CDec)).*cos(List1(i,CRA)-List2(:,CRA)));
      % how many objects has distances smaller then Teta
      Dist = sort(Dist);
      ObjNum = bin_sear(Dist,Teta);
      Cf(j,2) = Cf(j,2) + 2.*ObjNum./(N2.*(1-cos(Teta)));
   end
   Cf(j,1) = Teta;
   Cf(j,2) = Cf(j,2)./N1;
end




