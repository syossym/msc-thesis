%*************************************************************************
% QEBS 1.0: QCLs Electronic Band Structure
% The program is used to determine the electronic band structure of QCLs
% @Author: Le Quang Khai
% Email: ronaldokhai@yahoo.com
% Room 318 Woncheon Hall
% Electrical and Computer Engineering Dept.
% Ajou University
%*************************************************************************
clear

%L=23.3e-9;
L=76.5e-9
%L=48.7e-9

n=1000;
%x=25.4e-9;
x=0;
dx=(L-x)/(n-1);

nm=1e-9;
q=1.6e-19;
h=6.625e-34/(2.0*pi);
EigFile = fopen('EigValues.txt','w');

%********** Position dependent effective mass without non-parabolicity *****
a=createMatrixPosition(n,x,L,dx);
[v,e]=eig(a);
eeV=e./q;

LD=76.5e-9;
dy=LD/(n-1);
x=0:dx:L;

y=0;
for i=1:n
    Po(i)=V(y)/q;
    fprintf(EigFile,'%f \n',eeV(i,i));
    y=y+dy;
end
y=0:dy:LD;
plot(y,Po,'r');
hold;

for i=1:n
    eigdiagonal(i)=eeV(i,i);
end
[eigva,eigmin]=min(eigdiagonal)
s2=eigmin+1;
s3=eigmin+2;


v1=10*v(:,eigmin).^2+eigva;
v2=10*v(:,s2).^2+eigdiagonal(s2);
v3=10*v(:,s3).^2+eigdiagonal(s3);
v4=10*v(:,s3+1).^2+eigdiagonal(s3+1);
v5=10*v(:,s3+2).^2+eigdiagonal(s3+2);
v6=10*v(:,s3+3).^2+eigdiagonal(s3+3);
v7=10*v(:,s3+4).^2+eigdiagonal(s3+4);
v8=10*v(:,s3+5).^2+eigdiagonal(s3+5);
v9=10*v(:,s3+6).^2+eigdiagonal(s3+6);
v10=10*v(:,s3+7).^2+eigdiagonal(s3+7);
v11=10*v(:,s3+8).^2+eigdiagonal(s3+8);
v12=10*v(:,s3+9).^2+eigdiagonal(s3+9);
v13=10*v(:,s3+10).^2+eigdiagonal(s3+10);
v14=10*v(:,s3+11).^2+eigdiagonal(s3+11);

% **********************************
plot(x,v1);
plot(x,v2);
plot(x,v3);
plot(x,v4);
plot(x,v5);
plot(x,v6);
plot(x,v7);
plot(x,v8);
plot(x,v9);
plot(x,v10);
plot(x,v11);
plot(x,v12);
plot(x,v13);
plot(x,v14);
%grid on;
