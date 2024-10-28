%***************************************************************************************
% FEMQWEBS @Finite element method for quantum well electronic band structure calculation
% | ------------ |
% | Description: |
% | ------------ |
% The tool box provides the procedure to calculate all band edge energies 
% and corresponding wavefunctions in sinlge quantum square well using
% Finite Element Method.
% | ------------ |
% | Author:      |
% | ------------ |
% Full name: Le Quang Khai
% Email: ronaldokhai@yahoo.com or lqkhai@ajou.ac.kr
%***************************************************************************************
function [kk,ff]=feasmbl2(kk,ff,k,f,index)
% Purpose:
% Assembly of element matrices and vectors into the system matrix and vector
% Variable description
% kk: system matrix
% ff: system vector
% k: element matrix
% f: element vector
% index: dof vector associated with an element
edof=length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
        ff(ii,jj)=ff(ii,jj)+f(i,j);
    end
end