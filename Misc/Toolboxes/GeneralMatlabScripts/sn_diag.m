function [S,N]=sn_diag(Obs_Vec,Type,N_Bins);
%--------------------------------------------------------------------
% sn_diag function     Sum vs. Number diagram. Total population up to
%                    number vs. number.
% input  : - vector of observations.
%          - type of sort, 'a' - ascending.
%                          'd' - descending. (default).
%          - number of elements in the vector of numbers
%            (define the number axis).
%            default is [min:range/100:max].
% output : - vector of numbers (define the number axis).
%          - total number of objects with "observable value" up to
%            number defined by Sum_Vec.
%          + if number of output argument is smaller then 2 then
%            S-N diagram is plotted.
%    By  Eran O. Ofek           January 1998
%--------------------------------------------------------------------
if nargin==1,
   N_Bins = 100;
   Type = 'd';
   Sum_Vec = [max(Obs_Vec):-(max(Obs_Vec)-min(Obs_Vec))./N_Bins:min(Obs_Vec)]';
elseif nargin==2,
   N_Bins = 100;
   if Type=='d',
      Sum_Vec = [max(Obs_Vec):-(max(Obs_Vec)-min(Obs_Vec))./N_Bins:min(Obs_Vec)]';
   elseif Type=='a',
      Sum_Vec = [min(Obs_Vec):(max(Obs_Vec)-min(Obs_Vec))./N_Bins:max(Obs_Vec)]';
   else
      error('unknon Type');
   end
elseif nargin==3,
   % define Sum_Vec
   if Type=='d',
      Sum_Vec = [max(Obs_Vec):-(max(Obs_Vec)-min(Obs_Vec))./N_Bins:min(Obs_Vec)]';
   elseif Type=='a',
      Sum_Vec = [min(Obs_Vec):(max(Obs_Vec)-min(Obs_Vec))./N_Bins:max(Obs_Vec)]';
   else
      error('unknon Type');
   end
else
   error('number of argument must be 1,2 or 3');
end

Obs_Vec = sort(Obs_Vec);

S         = Sum_Vec;
N_Sum_Vec = length(Sum_Vec);
N         = zeros(N_Sum_Vec,1);

for Ind=1:1:N_Sum_Vec,
   N(Ind) = bin_sear(Obs_Vec,Sum_Vec(Ind));
end

if Type=='d',
   N = length(Obs_Vec) - N + 1;
elseif Type=='a',
   % do nothing
else
   error('unknon Type');
end


if (nargout==0 | nargout==1),
   plot(S,N);
end
