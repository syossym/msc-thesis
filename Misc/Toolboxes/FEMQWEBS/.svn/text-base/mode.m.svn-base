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
function m=mode(x)
nm=1e-9;
W=20e-9;
nmW=nm/W;
    if(x>=0.0&x<=10*nmW)
		m=1;
        return;
    end
	if(x>10*nmW&x<30*nmW)
        m=0;
		return;
    end
	if(x>=30*nmW&x<=40*nmW)
		m=1;
        return;
    end
    m=1;
    return ;
    