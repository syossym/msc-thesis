%***************************************************************************************
% FEMQWEBS @Finite element method for quantum well electronic band structure calculation
% | ------------ |
% | Description: |
% | ------------ |
% The tool box provides the procedure to calculate all band edge energies 
% and corresponding wavefunctions in sinlge quantum square well using
% Finite Element Method.
% | ------------ |
% | Author:      |
% | ------------ |
% Full name: Le Quang Khai
% Email: ronaldokhai@yahoo.com or lqkhai@ajou.ac.kr
%***************************************************************************************
clear
%-----------------------------------------------------------
% Input data for control parameters
%-----------------------------------------------------------
nel=1000;  % number of elements
nnel=2; % number of nodes per element
ndof=1; % number of dofs per node
nnode=1001;    % total number of nodes in system
sdof=nnode*ndof; % total system dofs
%-----------------------------------------------------------
% Input data for nodal coordinate values
%-----------------------------------------------------------
t=0;
L=40e-9     % Domain region
W=20e-9;    % Width of the quantum well
dx=L/(nnode-1);
for i=1:nnode
gcoord(i)=t;
t=t+dx/W;
%gcoord(i)=gcoord(i)/W;
end
    
%-----------------------------------------------------------
% Input data for nodal connectivity for each element
%-----------------------------------------------------------
% nodes(1,1)=1;nodes(1,2)=2;nodes(2,1)=2;nodes(2,2)=3;
% nodes(3,1)=3;nodes(3,2)=4;nodes(4,1)=4;nodes(4,2)=5;
% nodes(5,1)=5;nodes(5,2)=6;
for j=1:2
    for i=1:nel
        if(j==2)
            nodes(i,j)=i+1;
        else
            nodes(i,j)=i;
        end
    end
end
%-----------------------------------------------------------
% Input data for coefficients of the ODE
%-----------------------------------------------------------
% acoef=1;
% bcoef=-3;
% ccoef=2;
%-----------------------------------------------------------
% Input data for boundary conditions
%-----------------------------------------------------------
bcdof(1)=1; % first node is constrained
bcval(1)=0; % whose described value is 0
bcdof(2)=11; % sixth node is contraines
bcval(2)=0; % whose described value is 0
%-----------------------------------------------------------
% Initialization of matrices and vectors
%-----------------------------------------------------------
% for i=1:sdof
%     ff(i)=0;    % Initialization of system force vector
% end
%ff=zeros(sdof,sdof);
for i=1:sdof
    for j=1:sdof
        kk(i,j)=0;  % Initialization of system matrix
        ff(i,j)=0;  % Initialization of system force vector
    end
end
for i=1:nnel*ndof
    index(i)=0;     % Initialization of index vector
end
%-----------------------------------------------------------
% Compuation of element matrices and vectors and their assembly
%-----------------------------------------------------------
for i=1:nel                         % loop for the total number of elements
    nL=nodes(i,1);nR=nodes(i,2);    % extract nodes for ith element
    xL=gcoord(nL);xR=gcoord(nR);    % extract nodal coord values
    eleng=xR-xL;
    index=feeldof1(i,nnel,ndof);    % extraxt system dofs associated
    [k,p,q]=feMatrixK(xL,xR);             % compute element matrix
    f=feMatrixQ(xL,xR);             % compute element vector
    [kk,ff]=feasmbl2(kk,ff,k,f,index);  % assemble element matrices and vectors
end                                 % end of loop for total elements
%-----------------------------------------------------------
% Apply boundary conditions
%-----------------------------------------------------------
%[kk,ff]=feaplybc(kk,ff,bcdof,bcval);
%-----------------------------------------------------------
% Solve the matrix equation
%-----------------------------------------------------------
fem_eigma=ff^-1*kk;
[fem_eigvec,fem_eigva]=eig(fem_eigma);

h_=6.625e-34/(2*pi);
m0=9.1e-31;
mb=0.0901*m0;
Einf=h_^2*pi^2/(2*mb*W^2);
fem_eigva=(fem_eigva*Einf)/1.6e-19;

%***********************************************
%***** Plot the bound states and energy levels
y=0;
for i=1:nnode
    Po(i)=V(y)*Einf/1.6e-19;
    y=y+dx/W;
end
plot(gcoord,Po,'r');
hold;
%***********************************************
for i=1:nnode
    eigdiagonal(i)=fem_eigva(i,i);
end
[eigva,eigmin]=min(eigdiagonal)
s2=eigmin+1;
s3=eigmin+2;

v1=20*fem_eigvec(:,eigmin).^2+eigva;
v2=20*fem_eigvec(:,s2).^2+eigdiagonal(s2);
v3=20*fem_eigvec(:,s2+1).^2+eigdiagonal(s2+1);

%***********************************************
plot(gcoord,v1);
plot(gcoord,v2);
plot(gcoord,v3);


