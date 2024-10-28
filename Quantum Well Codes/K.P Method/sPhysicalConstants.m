%-------------------------------------------------------
%		
%   Name  : QW_K_P_PhysicalConstants 
%   Type  : Script
%	
%   Description : 										   
%	This file containts the common physical constants  
%	used by all the other scripts and functions in 
%   this project. All constants are given in CI units.
%
%   Usage : - 
%         
%
%   Tested : Matlab 7.5.0
%	    By : Yossi Michaeli, May 2009
%	
%-------------------------------------------------------

global s_PhysConsts;

s_PhysConsts.h         = 6.626176583e-34;    % Planck's constant [J sec]
s_PhysConsts.hbar      = 1.054588757e-34;    % Planck's constant over 2*pi [J sec]
s_PhysConsts.e_0       = 1.602189246e-19;    % electronic charge [C]
s_PhysConsts.m_0       = 9.109534e-31;       % electron rest mass [kg]
s_PhysConsts.c_0       = 2.99792461e+8;      % speed of light [m sec^-1]
s_PhysConsts.epsilon_0 = 8.854187827e-12;    % permittivity of free space [m^-3 kg^-1 sec^4 A^2]
s_PhysConsts.k_b       = 1.380662e-23;       % Boltzmann constant [m^2 kg s^-2 K^-1]
s_PhysConsts.mu_b      = 9.274e-24;          % Bohr magneton [J Tesla^-1]