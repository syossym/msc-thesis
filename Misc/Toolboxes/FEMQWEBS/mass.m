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
function [m]=mass(x)

nm=1e-9;
m0=9.1e-31;
    mw=0.0665*m0;
    mb=0.0901*m0;
    if(mode(x)==1)
        m=mb/mb;
        return
    else
        m=mw/mb;
        return
    end