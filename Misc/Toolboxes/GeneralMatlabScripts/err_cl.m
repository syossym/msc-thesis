function [Err,ErrErr]=err_cl(Data,Conf);
%--------------------------------------------------------------------
% err_cl function    Calculate right/left confidence boundary to a
%                  given numerical distribution.
% input  : - Data vector.
%          - column vector of confidence level,
%            default is [0.6827; 0.9545; 0.9973]
%            error is calculate to each one of the levels.
% output : - matrix of confidence interval (CI), the first column for
%            left CI and the second column to right CI.
%            Each line for each confidence level.
%          - estimate the standard error in the errors for
%            each confidence level.
% Tested : Matlab 5.3
%     By : Eran O. Ofek        November 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
if (nargin==1),
   Conf = [0.6827;0.9545;0.9973];
elseif (nargin==2),
   % do nothing
else
   error('Illigal number of input arguments');
end

N        = length(Data);
SortData = sort(Data);
Mean     = mean(Data);
Median   = median(Data);

BoundFrac = 0.5.*(1 - Conf);
NLower    = BoundFrac.*N;
NUpper    = (1-BoundFrac).*N;

% take mean value between floor an ceil
for I=1:1:length(Conf),
   if (floor(NLower(I))<1),
      LowerConfPos(I) = SortData(1);
   else
      LowerConfPos(I) = 0.5.*(SortData(floor(NLower(I))) + SortData(ceil(NLower(I))));
   end

   if (ceil(NUpper(I))>N),
      LowerConfPos(I) = SortData(N);
   else
      UpperConfPos(I) = 0.5.*(SortData(floor(NUpper(I))) + SortData(ceil(NUpper(I))));
   end
end
Err = [LowerConfPos', UpperConfPos'];

if (nargout==2),
   % poisson noise
   N_Err = sqrt(NLower);

   for I=1:1:length(Conf),
      if ((ceil(NLower(I)+N_Err(I))>N) | (floor(NLower(I)+N_Err(I))<1)),
         LowerConfPosEU = SortData(2);
      else
         LowerConfPosEU(I) = 0.5.*(SortData(floor(NLower(I)+N_Err(I))) + SortData(ceil(NLower(I)+N_Err(I))));
      end
   end

   for I=1:1:length(Conf),
      if (floor(NLower(I)-N_Err(I))<1),
         LowerConfPosEL(I) = SortData(1);
      else
         LowerConfPosEL(I) = 0.5.*(SortData(floor(NLower(I)-N_Err(I))) + SortData(ceil(NLower(I)-N_Err(I))));
      end
   end
   ErrErr = [LowerConfPosEU - LowerConfPosEL]';
end
