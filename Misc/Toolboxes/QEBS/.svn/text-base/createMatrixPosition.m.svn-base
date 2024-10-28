%*************************************************************************
% QEBS 1.0: QCLs Electronic Band Structure
% The program is used to determine the electronic band structure of QCLs
% @Author: Le Quang Khai
% Email: ronaldokhai@yahoo.com
% Room 318 Woncheon Hall
% Electrical and Computer Engineering Dept.
% Ajou University
%*************************************************************************
function [a] = createMatrixPosition(n,x,L,dx)
%L=23.3e-9
%     L=76.5e-9;
%      x=0.0;
%      dx=L/(n-1);
h=6.625e-34/(2.0*pi);

alpha_w=-h*h/(dx*dx*(mass(x)+mass(x-dx)));
alpha_e=-h*h/(dx*dx*(mass(x+dx)+mass(x)));
alpha_x=-alpha_w-alpha_e;

a(1,1)=alpha_x+V(x);
a(1,2)=alpha_e;

max=n;
for i=1:n
    a(3,i)=0.0;
end
for i=1:n-2
    a(max,i)=0.0;
end
for i=2:n-1
    x=x+dx;
    for j=2:n-1
        if(j~=i) a(i,j+1)=0.0; end
        if(j==i)
            alpha_w=-h*h/(dx*dx*(mass(x)+mass(x-dx)));
            alpha_e=-h*h/(dx*dx*(mass(x+dx)+mass(x)));
            alpha_x=-alpha_w-alpha_e;
            a(i,j)=alpha_x+V(x);
            a(i,j-1)=alpha_w;
            a(i,j+1)=alpha_e;
        end
    end
end
alpha_w=-h*h/(dx*dx*(mass(x)+mass(x-dx)));
alpha_e=-h*h/(dx*dx*(mass(x+dx)+mass(x)));
alpha_x=-alpha_w-alpha_e;
a(max,max-1)=alpha_w;
a(max,max)=alpha_x+V(x);

