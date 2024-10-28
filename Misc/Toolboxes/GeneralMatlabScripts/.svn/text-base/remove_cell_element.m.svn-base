function NewCell=remove_cell_element(Cell,El2Remove);
%-------------------------------------------------------------------------
% remove_cell_element function   Remove a list of indices from a cell
%                              vector.
% Input  : - Cell vector.
%          - Indices of elements to remove from cell vector.
% Output : - New cell vector..
% Tested : Matlab 7.0
%     By : Eran O. Ofek             June 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------

N2rem   = length(El2Remove);
Ncell   = length(Cell);

NewInd  = setdiff([1:1:Ncell].',El2Remove);
NewCell = Cell(NewInd);

