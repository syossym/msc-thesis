function [ErrExp,ErrVar]=symerror(Expression,varargin);
%-------------------------------------------------------------------
% symerror function     Given a symbolic expression and the
%                     variables of the expression, calculate
%                     the symbolic error function of the
%                     expression, with respect to the variables.
%                     The new expression contains an error
%                     variable named D_"original_var" for each
%                     variable respectively.
%                     Use: char(vectorize(ErrExp)) to convert
%                     the symbolic expression to a vectorized
%                     string that can be evaluated using eval.
%                     The function also returns a cell array of
%                     the new error variables
%                     (e.g.,  D_"original_var").
% Input  : - Symbolic expression.
%          * Arbitrary number of variables to differentiate by.
% Output : - Symbolic error expression.
%          - Cell array of new error variables.
% Tested : Matlab 6.5
%     By : Eran O. Ofek            August 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [ErrExp,ErrVar]=symerror('sqrt(x^2+y^2)','x','y')
%-------------------------------------------------------------------

Nvar = length(varargin);
for I=1:1:Nvar,
   Var         = varargin{I};
   eval(sprintf('syms D_%s',char(Var)));
   ErrVar{I}   = sprintf('D_%s',char(Var));
   SubDiff{I}  = (diff(Expression,Var)*sprintf('D_%s',char(Var)))^2;
end

syms ErrExp;
ErrExp = 0;
for I=1:1:Nvar,
   ErrExp = ErrExp + SubDiff{I};
end
ErrExp = sqrt(ErrExp); 
