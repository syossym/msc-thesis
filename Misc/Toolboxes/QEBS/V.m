%*************************************************************************
% QEBS 1.0: QCLs Electronic Band Structure
% The program is used to determine the electronic band structure of QCLs
% @Author: Le Quang Khai
% Email: ronaldokhai@yahoo.com
% Room 318 Woncheon Hall
% Electrical and Computer Engineering Dept.
% Ajou University
%*************************************************************************
function [v]=V(x)
V0=0.295*1.6e-19;
nm=1e-9;
F=-48e+5;
q=1.6e-19;
if(mode(x)==1)
    v=V0+q*F*x;
    return
else
    v=q*F*x;
    return
end
