function V=mat2vec(M);
%--------------------------------------------------------------------
% mat2vec function     convert matrix to vector.
% Input  : - Matrix.
% Output : - Vector.
% Tested : Matlab 5.3
%     By : Eran O. Ofek             May 2000   
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%--------------------------------------------------------------------

[SizeI,SizeJ] = size(M);

V = zeros(SizeI.*SizeJ,1);

for K=1:1:SizeJ,
   V((K-1).*SizeI+1:K.*SizeI) = M(:,K);
end

