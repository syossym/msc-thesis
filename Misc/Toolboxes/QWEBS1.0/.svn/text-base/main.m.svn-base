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
clear

%step=99
%n=2*step+3+1;
n=501
z=0;
q=1.6e-19;
L=20e-9;
%nnode=step+1;
%dz=1/nnode;
dz=L/(n-1);
a=create2x2UMatrix(n,L,z,dz);

h_=6.625e-34/(2*pi);
m0=9.1e-31;
Einf=h_^2/(m0*L^2);

[v,e]=eig(a);
eeV=e./q;

y=0;
for i=1:n
    Po(i)=V(y,0)/q;
    y=y+dz;
end
    yy=0:dz:L;
figure;plot(yy,Po,'r');
hold

for i=1:2*n-2
    eigdiagonal(i)=eeV(i,i);
end
% for i=1:n-2
%     eigdiagonal(i)=eeV(i,i);
% end
[eigva,eigmin]=maxneg(eigdiagonal,2*n-2)
s2=eigmin+1;
s3=eigmin+2;


v1=5*v(:,eigmin).^2+eigva;
v2=5*v(:,s2).^2+eigdiagonal(s2);
v3=5*v(:,s3).^2+eigdiagonal(s3);
v4=5*v(:,s3+1).^2+eigdiagonal(s3+1);
% v5=20*v(:,s3+2).^2+eigdiagonal(s3+2);
% v6=20*v(:,s3+3).^2+eigdiagonal(s3+3);
% **********************************
y=1:2*n-2;
figure;plot(y,v1);
hold;
plot(y,v2);
plot(y,v3);
plot(y,v4);
% plot(y,v5);
% plot(y,v6);