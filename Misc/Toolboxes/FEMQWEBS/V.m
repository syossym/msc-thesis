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
function [v]=V(x)
% Potential without Hartree effect
V0=0.276*1.6e-19;
nm=1e-9;
h_=6.625e-34/(2*pi);
m0=9.1e-31;
mb=0.0901*m0;
L=20e-9;
%F=-50e+5;
F=0.0;
q=1.6e-19;
Einf=h_^2*pi^2/(2*mb*L^2);
if(mode(x)==1)
    v=(V0+q*F*x)/Einf;
    return  
else
    v=q*F*x/Einf;
    return
end