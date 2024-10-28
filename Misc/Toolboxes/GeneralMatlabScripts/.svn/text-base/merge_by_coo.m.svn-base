function MerCat=merge_by_coo(Cat1,Cat2,Threshold,DistMethod,Cat1Col,Cat2Col);
%------------------------------------------------------------------
% merge_by_coo function     merge two catalogs by coordinates.
%                         assume coordinates are in radians.
%                         reject entries that do not have a match.
% Input  : - First catalog.
%          - Second catalog.
%          - merging threshold in arcsec/pixels.
%            if the distnace method is 's' then arcsec, else pixels.
%          - Distance calculating method:
%            's' - spherical triangle (coordinated in radians) - default.
%            'p' - pitagorian.
%          - First catalog coordinates columns, default is [1 2].
%          - Second catalog coordinates columns, default is [1 2].
% Output : - Merged catalog. [Cat1, Cat2]
% Tested : Matlab 5.3
%     By : Eran O. Ofek          March 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------
RADIAN = 57.29577951308232;
if (nargin==3),
   DistMethod='s';
   Cat1Col = [1 2];
   Cat2Col = [1 2];
elseif (nargin==4),
   Cat1Col = [1 2];
   Cat2Col = [1 2];
elseif (nargin==5),
   Cat2Col = [1 2];
elseif (nargin==6),
   % do nothing
else
   error('illegal number of input arguments');
end

% convert threshold from arcsec to radians
if (DistMethod=='s'),
   Threshold = Threshold./(3600*RADIAN);
else
   Threshold = Threshold;
end

[N1, C1] = size(Cat1);
[N2, C2] = size(Cat2);

NewCat = zeros(max([N1,N2]),C1+C2);

J = 0;
for I=1:1:N1,
   RA1  = Cat1(I,Cat1Col(1));
   Dec1 = Cat1(I,Cat1Col(2));
   RA2  = Cat2(:,Cat2Col(1));
   Dec2 = Cat2(:,Cat2Col(2));

   if (DistMethod=='s'),
      Dist = sphere_dist(RA1,Dec1,RA2,Dec2);
   elseif (DistMethod=='p'),
      Dist = sqrt((RA1-RA2).^2 + (Dec1-Dec2).^2);
   else
      error('unknown distance method');
   end
   IfD = find(Dist<=Threshold);
   if (isempty(IfD)),
      % counterpart not found
   else
      % counterpart found, take closest
      [MinDist, MinDistI] = min(Dist);
      J = J + 1;
      NewCat(J,:) = [Cat1(I,:), Cat2(MinDistI,:)];
   end
end

MerCat = NewCat(1:J,:);
