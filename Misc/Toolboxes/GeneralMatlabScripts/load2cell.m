function [CellMat]=load2cell(FileName,varargin);
%----------------------------------------------------------------------------
% load2cell function     Read table into cell array.
%                      This is useful whenever there is text
%                      in the table.  
%                      Numbers are stored as numbers and empty 'cells' are
%                      filled with NaNs.
% Input  : - File name to load.
%          * Pairs of parameters ...,keyword,value,...
%            where keywords are:
%            'delim'   - delimiter, default is 'spaces'
%                        Note there is a diiference between ' '
%                        and 'spaces'.
%            'skipN'   - Number of header line to skip, default is 0.
%            'skipC'   - skip lines starting with a given character,
%                        default is '%'
% Output : - Cell array.
% Tested : Matlab 5.3           
%     By : Eran O. Ofek            August 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------

Narg = length(varargin);
if (Narg./2~=floor(Narg./2)),
   error('Illegal number of input arguments');
end

Delimiter   = 'spaces';
SkipLineN   = 0;
SkipLineC   = '%';

for I=1:2:(Narg-1),
   switch varargin{I}
    case 'delim'
       Delimiter   = varargin{I+1};
    case 'skipN'
       SkipLineN   = varargin{I+1};
    case 'skipC'
       SkipLineC   = varargin{I+1};
    otherwise
       error('Unkonwn keyword parameter');
   end
end

FID = fopen(FileName,'r');

for I=1:1:SkipLineN,
   Junk = fgetl(FID);
end

LenSkipLineChar = length(SkipLineC);

CellMat = cell(0,0);
K = 0;
while (feof(FID)==0),
   Line = fgetl(FID);

   if (strcmp(Line(1:LenSkipLineChar),SkipLineC)==1),
      % skip line
   else
      % read line
      K = K + 1;

      switch Delimiter
       case 'spaces'
          while (~isempty(findstr(Line,'  '))),
             Line = strrep(Line,'  ',' ');
          end
       otherwise
          % do nothing
      end
   
      switch Delimiter
       case 'spaces'
          Ind  = [0, findstr(Line,' '), length(Line)+1];
       otherwise
          Ind  = [0, findstr(Line,Delimiter), length(Line)+1];
      end 
      N    = length(Ind);
   
   
      Cell_Col = 0;
      for I=1:1:N-1,
         Cell_Col = Cell_Col + 1;
         String = Line(Ind(I)+1:Ind(I+1)-1);
         if (isempty(String)==1),
            CellMat{K,Cell_Col} = NaN;
         else
            if (isempty(str2num(String))==1),
               CellMat{K,Cell_Col} = String;
            else
               CellMat{K,Cell_Col} = str2num(String);
            end
         end
      end
   end
end   
