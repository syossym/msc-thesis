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
function [maxneg_eigva,maxneg_eigindex]=maxneg(a,n)

for i=1:n
    if(a(i)<0)
        init_t=a(i);        
        index=i;
        break
    end
end
for(i=index:n)
    if(a(i)<0&&abs(init_t)>abs(a(i)))
        init_t=a(i);
        index_max=i;
    end
end
maxneg_eigva=init_t;
maxneg_eigindex=index_max;
return