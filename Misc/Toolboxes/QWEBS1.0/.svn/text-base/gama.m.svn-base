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
function [g]=gama(x,mo)
g1GaAs=6.85;g2GaAs=2.1;g3GaAs=2.9;
%g1AlGaAs=6.85;g2AlGaAs=2.1;g3AlGaAs=2.9;
g1AlGaAs=5.29;g2AlGaAs=1.72;g3AlGaAs=2.46;

    if(mode(x)==0)
        if(mo==1)
            g=g1GaAs;
            return
        elseif(mo==2)
            g=g2GaAs;
            return
        else
            g=g3GaAs;
            return
        end    
    else
        if(mo==1)
            g=g1AlGaAs;
            return
        elseif(mo==2)
            g=g2AlGaAs;
            return
        else
            g=g3AlGaAs;
            return
        end    
    end