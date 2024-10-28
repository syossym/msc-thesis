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
function [kk,ff]=feaplybc(kk,ff,bcdof,bcval)
% Purpose:
% Apply constraints to matrix equation [kk]x=lamda[ff]x
n=length(bcdof);
sdof=size(kk);
for i=1:n
    c=bcdof(i);
    for j=1:sdof
        kk(c,j)=0;
    end
    kk(c,c)=1;
    ff(c)=bcval(i);
end