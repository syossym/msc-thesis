%*************************************************************************
% QEBS 1.0: QCLs Electronic Band Structure
% The program is used to determine the electronic band structure of QCLs
% @Author: Le Quang Khai
% Email: ronaldokhai@yahoo.com
% Room 318 Woncheon Hall
% Electrical and Computer Engineering Dept.
% Ajou University
%*************************************************************************
function m=mode(x)
nm=1e-9;

% The period of QC laser including two injector and one active region.

if(x>=0.0&x<=3.4*nm)
    m=1;
    return;
end

if(x>3.4*nm&x<6.6*nm)
    m=0;
    return;
end
if(x>=6.6*nm&x<=8.6*nm)
    m=1;
    return;
end
if(x>8.6*nm&x<11.4*nm)
    m=0;
    return;
end
if(x>=11.4*nm&x<=13.7*nm)
    m=1;
    return;
end
if(x>13.7*nm&x<16.0*nm)
    m=0;
    return;
end
if(x>=16.0*nm&x<=18.5*nm)
    m=1;
    return;
end
if(x>18.5*nm&x<20.8*nm)
    m=0;
    return;
end
if(x>=20.8*nm&x<=23.3*nm)
    m=1;
    return;
end
if(x>23.3*nm&x<25.4*nm)
    m=0;
    return;
end
if(x>=25.4*nm&x<=31.2*nm)
    m=1;
    return;
end
if(x>31.2*nm&x<32.7*nm)
    m=0;
    return;
end
if(x>=32.7*nm&x<=34.7*nm)
    m=1;
    return;
end
if(x>34.7*nm&x<39.6*nm)
    m=0;
    return;
end
if(x>=39.6*nm&x<=41.3*nm)
    m=1;
    return;
end
if(x>41.3*nm&x<45.3*nm)
    m=0;
    return;
end
if(x>=45.3*nm&x<=48.7*nm)
    m=1;
    return;
end
if(x>48.7*nm&x<51.9*nm)
    m=0;
    return;
end
if(x>=51.9*nm&x<=53.9*nm)
    m=1;
    return;
end
if(x>53.9*nm&x<56.7*nm)
    m=0;
    return;
end
if(x>=56.7*nm&x<=59.0*nm)
    m=1;
    return;
end
if(x>59.0*nm&x<61.3*nm)
    m=0;
    return;
end
if(x>=61.3*nm&x<=63.8*nm)
    m=1;
    return;
end
if(x>63.8*nm&x<66.1*nm)
    m=0;
    return;
end
if(x>=66.1*nm&x<=68.6*nm)
    m=1;
    return;
end
if(x>68.6*nm&x<70.7*nm)
    m=0;
    return;
end
if(x>=70.7*nm&x<=76.5*nm)
    m=1;
    return;
end

m=1;
return;