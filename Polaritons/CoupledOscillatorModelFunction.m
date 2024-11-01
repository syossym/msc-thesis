function [D_abs, D_real, D_imag, V_vec] = CoupledOscillatorModelFunction(params, E)
%
% This function is the target function for Couple Oscillator model 
% fitting procedure.
%
%   Input:  'params' - the model fitted parameters list.
%           'E' - the cavity mode energy [eV].
%
%   Output: 'F' - model fitting results - the Polariton branches.
%
% Tested: Matlab 7.6.0
% Created by: Yossi Michaeli, February 2010
% Edited by: -
%

global Consts;

%% Init the fitting parameters

% E_x1 = params(1);            % [eV]
% E_x2 = params(2);            % [eV]
% E_Omega_1 = params(3);       % [eV]
% E_Omega_2 = params(4);       % [eV]
% gamma_MC = params(5);        % [eV]
% gamma_x1 = params(6);        % [eV]
% gamma_x2 = params(7);        % [eV]

N = length(params(2:end))/3;
gamma_MC = params(1);
gamma_X = params(2:N+1);
E_Omega = params(N+2:2*N+1);
E_X  = params(2*N+2:3*N+1);

%% Polaritons branches calculation

D_real = [];
D_imag = [];
D_abs = [];
V_vec = cell(1,1+length(E_X));
for (dd=1:length(E))
    H_mat = zeros(1+length(E_Omega));
    H_mat(1,1) = E(dd)-1i*gamma_MC;
    H_mat = H_mat + diag([0, E_X-1i.*gamma_X]);
    H_mat(1,2:end) = E_Omega/2;
    H_mat(2:end,1) = E_Omega'/2;
    [V,D] = eig(H_mat);
    D_diag = diag(D);
    [temp1, I1] = sort(real(diag(D)));
    [temp2, I2] = sort(imag(diag(D)));
    [temp3, I3] = sort(abs(diag(D)));
    V_sorted = V(:,I3);
    imag_D_sorted = imag(D_diag(I3));
    real_D_sorted = real(D_diag(I3));
    %V_sroted = V_sorted.^2;
    D_real = [D_real, real_D_sorted];
    D_imag = [D_imag, imag_D_sorted];
    D_abs = [D_abs, temp3];
    for (vv=1:length(temp3))
        V_vec{vv} = [V_vec{vv}, V_sorted(:,vv)];
    end
end
% for (dd=1:length(E))
%     H_mat = [E(dd)+1i*gamma_MC, E_Omega_1/2, E_Omega_2/2 ;...
%              conj(E_Omega_1/2), E_x1+1i*gamma_x1, 0 ;...
%              conj(E_Omega_2/2), 0, E_x2+1i*gamma_x2];
%     [V,D] = eig(H_mat);
%     temp = sort(real(diag(D)));
%     F = [F, temp];
% end

