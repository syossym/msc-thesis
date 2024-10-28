function tf = strlexgt(a, b)
%STRLEXGT Lexicographic greater than.
%
%   STRLEXGT(A, B) returns 1 if A is lexicographically greater than B, and 0
%   otherwise.
%
%   This is a MATLAB version of the Perl `gt' operator.
%
%   See also GT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 11:15:17 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   tf = strlexcmp(a, b) > 0;
