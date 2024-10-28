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
function [k,p,q]=feMatrixK(xL,xR)
% Purpose 
% element matrix for [P] + pi^2*[Q] using liner element
eleng=xR-xL;
xav=(xR+xL)/2;
a=1/(mass(xav)*eleng);
b=pi^2*V(xav)*eleng/6;
k=[a+2*b -a+b; -a+b a+2*b];
p=[a -a;-a a];
q=[2*b b; b 2*b];