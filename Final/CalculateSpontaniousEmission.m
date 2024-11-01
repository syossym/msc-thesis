function QWParams = CalculateSpontaniousEmission(QWParams, Params)
%
% This function calculates the spontanious emission spectrum.
%
%   Input: 'QWParams' - structure containing calculated optical
%                       parameters for the simulated structure.
%          'Params' - simulaion parameters.
%
%   Output: 'QWParams' - structure containing calculated optical
%                       parameters for the simulated structure.
%
% Tested: Matlab 7.10.0
% Created by: Yossi Michaeli, September 2010
% Edited by: -
%

global Consts;

[n,m] = size(QWParams);

for (nn=1:n) 
    for (mm=1:m)
        QWParams{nn,mm}.Rsp_new.TE.HF = (1./(1-exp((Params.E_exc-(QWParams{nn,mm}.E_fc-QWParams{nn,mm}.E_fv))./(Consts.k_B.*Params.T)))).*...
                 QWParams{nn,mm}.alpha.TE.HF.*Params.E_exc.^2;
        QWParams{nn,mm}.Rsp_new.TM.HF = (1./(1-exp((Params.E_exc-(QWParams{nn,mm}.E_fc-QWParams{nn,mm}.E_fv))./(Consts.k_B.*Params.T)))).*...
                 QWParams{nn,mm}.alpha.TM.HF.*Params.E_exc.^2;     
        QWParams{nn,mm}.Rsp_new.HF = (2/3).*QWParams{nn,mm}.Rsp_new.TE.HF + (1/3).*QWParams{nn,mm}.Rsp_new.TM.HF;
        QWParams{nn,mm}.Rsp_new.TE.FCT = (1./(1-exp((Params.E_exc-(QWParams{nn,mm}.E_fc-QWParams{nn,mm}.E_fv))./(Consts.k_B.*Params.T)))).*...
                 QWParams{nn,mm}.alpha.TE.FCT.*Params.E_exc.^2;
        QWParams{nn,mm}.Rsp_new.TM.FCT = (1./(1-exp((Params.E_exc-(QWParams{nn,mm}.E_fc-QWParams{nn,mm}.E_fv))./(Consts.k_B.*Params.T)))).*...
                 QWParams{nn,mm}.alpha.TM.FCT.*Params.E_exc.^2;  
        QWParams{nn,mm}.Rsp_new.FCT = (2/3).*QWParams{nn,mm}.Rsp_new.TE.FCT + (1/3).*QWParams{nn,mm}.Rsp_new.TM.FCT;
    end
end