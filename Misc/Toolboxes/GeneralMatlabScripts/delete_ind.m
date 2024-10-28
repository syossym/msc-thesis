function NewMat=delete_ind(Mat,Ind,Dim);
%---------------------------------------------------------------------------
% insert_ind function    Delete a column/s or row/s from a specific
%                      position in a matrix.
% Input  : - Matrix.
%          - Indices of rows or columns to remove.
%          - Dimension: 1 - Delete rows; 2 - Delete columns, default is 1.
% Tested : Matlab 7.0
%     By : Eran O. Ofek             December 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: insert_ind.m
% Example: NewMat=insert_val(zeros(3,2),2,[1 1; 2 2],1)
%---------------------------------------------------------------------------
if (nargin==2),
   Dim = 1;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


Nind   = length(Ind);

Nsize  = size(Mat,Dim);
IndVec = [1:1:Nsize];
for Iind=1:1:Nind,
   IndVec(find(IndVec==Ind(Iind))) = NaN;   % set indices to delete to NaN
end
NewIndVec = IndVec(find(isnan(IndVec)==0));

switch Dim
 case 1
    NewMat = Mat(NewIndVec.',:);
 case 2
    NewMat = Mat(:,NewIndVec);
 otherwise
   error(sprintf('%d-dimension is unspported - use only 1/2-d',Dim));
end


