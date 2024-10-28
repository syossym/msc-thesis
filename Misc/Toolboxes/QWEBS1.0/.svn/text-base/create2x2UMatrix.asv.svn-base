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
function [ma] = create2x2UMatrix(n,L,z,dz)
h_=6.625e-34/(2.0*pi);
m0=9.1e-31;
aInGaSb=0.6077588e-9;
%kp=0.02;
kp=0;
%Initize a system matrix
for i=1:2*n-2
    for j=1:2*n-2
        ma(i,j)=0;
    end
end
%      for i=1:n-1
%          for j=1:n-1
%             ma(i,j)=0;
%          end
%      end
%Ec=Ev+Eg for an unstrained semiconductor
%********************************************************************
A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
%********************************************************************

a11=A1+V(z,0)+2*A2/(dz*dz);a12=-A2/(dz*dz);a1n=C;a1np1=-D/(2*dz);

ma(1,1)=a11;ma(1,2)=a12;ma(1,n)=a1n;ma(1,n+1)=a1np1;

for i=2:n-2
    z=z+dz;
    %********************************************************************
    A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
    A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
    B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
    B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
    C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
    D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
    %********************************************************************
    for j=2:n-2
        if(j==i)
            ma(i,j)=A1+2*A2/(dz*dz)+V(z,0);
            ma(i,j-1)=-A2/(dz*dz);
            ma(i,j+1)=-A2/(dz*dz);
        end
    end
    for j=n+1:2*n-3
        if(j==n+i-1)
            ma(i,j)=C;
            ma(i,j-1)=D/(2*dz);
            ma(i,j+1)=-D/(2*dz);
        end
    end
end
z=L;
%********************************************************************
A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
%********************************************************************
ma(n-1,n-2)=-A2/(dz*dz);ma(n-1,n-1)=A1+2*A2/(dz*dz)+V(z,0);
ma(n-1,2*n-3)=D/(2*dz);ma(n-1,2*n-2)=C;

z=0;
%********************************************************************
A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
%********************************************************************
ma(n,1)=C;
ma(n,2)=D/(2*dz);
ma(n,n)=B1+B2/(dz*dz)+V(z,0);
ma(n,n+1)=-B2/(dz*dz);

for i=n+1:2*n-3
    z=z+dz;
    %********************************************************************
    A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
    A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
    B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
    B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
    C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
    D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
    %********************************************************************
    for j=2:n-2
        if(j==i-n+1)
            ma(i,j)=C;
            ma(i,j-1)=-D/(2*dz);
            ma(i,j+1)=D/(2*dz);
        end
    end
    for j=n+1:2*n-3
        if(j==i)
            ma(i,j)=B1+2*B2/(dz*dz)+V(z,0);
            ma(i,j-1)=-B2/(dz*dz);
            ma(i,j+1)=-B2/(dz*dz);
        end
    end
end
z=L;
%********************************************************************
A1=0.5*(gama(z,1)+gama(z,2))*kp^2*h_^2/m0;
A2=-0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
B1=0.5*(gama(z,1)-gama(z,2))*kp^2*h_^2/m0;
B2=-0.5*(gama(z,1)+2*gama(z,2))*h_^2/m0;
C=sqrt(3)*(gama(z,2)+gama(z,3))*kp^2*h_^2/(4*m0);
D=sqrt(3)*gama(z,3)*kp*h_^2/m0;
%********************************************************************
ma(2*n-2,n-1)=C;ma(2*n-2,n-2)=-D/(2*dz);
ma(2*n-2,2*n-3)=-B2/(dz*dz);ma(2*n-2,2*n-2)=B1+2*B2/(dz*dz)+V(z,0);


%      A2=0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
%
%     %********************************************************************
%
%      a11=V(z,0)-2*A2/(dz*dz);a12=A2/(dz*dz);
%
%      ma(1,1)=a11;ma(1,2)=a12;
% %      ma(1,n)=a1n;ma(1,n+1)=a1np1;
%
%      for i=2:n-1
%          z=z+dz;
%          %********************************************************************
%           A2=0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
%          for j=2:n-1
%              if(j==i)
%                  ma(i,j)=-2*A2/(dz*dz)+V(z,0);
%                  ma(i,j-1)=A2/(dz*dz);
%                  ma(i,j+1)=A2/(dz*dz);
%              end
%          end
% %          for j=n+1:2*n-3
% %              if(j==n+i-1)
% %                  ma(i,j)=C;
% %                  ma(i,j-1)=D/(2*dz);
% %                  ma(i,j+1)=-D/(2*dz);
% %              end
% %          end
%      end
%      %********************************************************************
%     A2=0.5*(gama(z,1)-2*gama(z,2))*h_^2/m0;
%     %********************************************************************
%      ma(n,n-1)=A2/(dz*dz);ma(n,n)=-2*A2/(dz*dz)+V(z,0);
