%1D-simulation with AQUILA
%modulation doped high-mobility sample

initaquila                         %initialize AQUILA
aquila_control.mode=1;             %1D-Simulation

%structure definition starting at surface
add_mbox(100,20,0,0);              %100A GaAs cap
add_mbox(1490,100,0.3,0);          %1500A AlGaAs
add_mbox(3,1,0.3,0);               %3A AlGaAs around doping with enhanced resolution
add_mbox(4,1,0.3,-2.09e19*0.94);   %4A delta-doping in AlGaAs
add_mbox(3,1,0.3,0);               %3A AlGaAs with enhanced resolution
add_mbox(800,20,0.3,0);            %800A AlGaAs spacer
add_mbox(12000,100,0.0,0);         %more than 10000A GaAs, transition to bulk

add_qbox([2350 3500],10,3,GE+HH+LH);     %the quantum region, where the 2DEG resides

add_pbox([2350 3500], CB+VB);          %zoom in the channel region
add_pbox([1590 1600], CB+VB);          %monitor the doping layer for parallel channel
add_pbox([0 12000], CB+VB);   

add_boundary(LEFT,POTENTIAL,0);    %surface, no gate
add_boundary(RIGHT,FIELD,0);       %transition to bulk

runstructure;                      %DO IT
