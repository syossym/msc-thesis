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
function [f]=feMatrixQ(xL,xR)
% Purpose:
% element vector for pi^2*[R] using liner element
eleng=xR-xL;    % element length
a=pi^2*eleng/6;
f=[2*a a; a 2*a];

