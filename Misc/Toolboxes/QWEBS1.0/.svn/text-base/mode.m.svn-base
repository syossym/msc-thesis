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
function m=mode(x)
nm=1e-9;
L=20*nm;
nmW=nm/L;
% The period of QC laser including two injector and one active region.

    po=[5,10,5];
    x1=po(1);x2=po(2);x3=po(3);
    if(x>=0.0&x<=x1*nm)
		m=1;  
        return;
    end
	if(x>x1*nm&x<(x1+x2)*nm)
		m=0;
       return;
    end
	if(x>=(x1+x2)*nm&x<=(x1+x2+x3)*nm)
		m=1;
        return;
    end    
    m=1;
    return;