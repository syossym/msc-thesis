Simple example about how to calculate the circulating power in Fabry Perot cavity
To run OSCAR, execute the script "Run_OSCAR.m"

New version, optimised for speed:
- New method to find the resonance length. The length of the cavity is adjusted to ensure the circulating field is added coherently. So it is no longer needed to scan the cavity power over one FSR
- Change the definition of 'Mat_propagation' to not have to use 'fftshift' procedure
- Precalculate the wavefront induced by the mirrors.
- Added a procedure to calculate the number of iteration as a function of the desired accuracy.