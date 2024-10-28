%*************************************************************************
% QWEBS 1.0: Quantum Well Electronic Band Structure 
% The program is used to determine the electronic valence band structure of single
% quantum well using 4x4 k.p method.
% Example: In this code, the valence band of GaAs/AlGaAs single quantum
% well is determined.
% @Author: Le Quang Khai
% Email: ronaldokhai@yahoo.com
% Room 318 Woncheon Hall
% Electrical and Computer Engineering Dept.
% Ajou University
%*************************************************************************
function [v]=V(x,type)
V0e=0.2506*1.6e-19;

V0h=0.1235*1.6e-19;

nm=1e-9;
m0=9.1e-31;
%F=-50e+5;
F=0;
q=1.6e-19;
h_=6.625e-34/(2*pi);
L=20*nm;
Einf=h_^2/(m0*L^2);
    if(type==1) % the particle is an electron
        if(mode(x)==1)
            v=V0e+q*F*x;
            return
        else
            v=q*F*x;
            return
        end
    else % the particle is a hole
        if(mode(x)==1)
            v=(-V0h+q*F*x);
            return
        else
            v=q*F*x;
            return
        end
    end
  
