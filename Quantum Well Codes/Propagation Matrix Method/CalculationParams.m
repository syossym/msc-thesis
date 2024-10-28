% This script file contains some parameters for the 
% dispersion calculation and band-edge calculation 
% scripts.

%% Common params

min_en_diff = 0.05*meV;    % helps to determine whether two solutions are the same
a_lim       = 1e5;         % max. for a
k_t         = 1e5;         % k_t for band-edge calculation
E_step      = 0.1*meV;     % energy step for band-edge calculation
E_acc       = meV*1e-9;    % energy resolution of the calculated subband energies

%% Dispersion calculation - search parameters

subband_en_search = 0.15*meV;  % energy search range
en_abort          = meV;       % abort search when E range > en_abort
k_t_max           = 1e9;       % max. calculated k_t
k_t_step          = 2e6;       % k_t step for E-k dispersion
num_cb_subbands   = 4;         % number of conduction band subbands taken into account
num_states        = nr-1;      % number of valence band subbands taken into account