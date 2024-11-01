function [Res,IndLine]=mgrep(File,LookStr,varargin);
%------------------------------------------------------------------
% mgrep function     grep-like utility for MATLAB.
%                  Search for substrings within text file.
% Input  : - File name.
%          - lookup string.
%          * Arbitrary number of pairs of arguments:
%            ...,keyword,value,...
%            Available keywords are:
%            'CaseS'   - case sensitive {'y' | 'n'}, default is 'y'.
% Output : - Cell array of all lines containing look up string.
%          - Vector of positions of lookup string within line.
% Tested : Matlab 7.0
%     By : Eran O. Ofek        Feb 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------

Narg = length(varargin);
if (Narg.*0.5~=floor(0.5.*Narg)),
   error('Illegal number of input arguments');
end
% default arguments
CaseS     = 'y';

for Iarg=1:2:Narg-1,
   switch varargin{Iarg}
    case 'CaseS'
       CaseS     = varargin{Iarg+1};
    otherwise
       error('Unknown keyword Option');
   end
end

FID = fopen(File,'r');
Counter  = 0;
Continue = 1;
while (Continue==1),
   if feof(FID), 
      Continue = 0;
   else
      Line = fgetl(FID);
      Ind  = strfind(Line,LookStr);
      if (isempty(Ind)==0),
 	 Counter           = Counter + 1;
 	 Res{Counter}      = Line;
         IndLine(Counter)  = Ind;
      end
   end
end
fclose(FID);
