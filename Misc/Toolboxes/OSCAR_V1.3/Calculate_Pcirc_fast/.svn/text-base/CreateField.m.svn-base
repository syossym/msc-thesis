clear all

global Grid
global Laser
global Mirror
global ITM
global ETM
global Length;
global Mat_propagation
global Field


%-----------------------------------------------------------------
fprintf('\n*******************************************************\n ')
fprintf('                   OSCAR FFT code                       \n ')

%--------------------- Parameters simulations --------------------

Grid.Num_point = 128;
Grid.Length = 0.30;

% Parameters of the input laser beam
Laser.lambda = 1064e-9;
Laser.beam_radius = 0.0200;
Laser.wavefront_radius = -2000;

Laser.power = 1; 

Laser.amplitude = sqrt(Laser.power*2/(pi*Laser.beam_radius^2)); % Amplitude for TEM00 input beam
%--------------------- Cavity parameters ---------------------------

Length_cav = 1000;

ITM.RofC = 2500;
ETM.RofC = 2500;

Mirror.Diam = 0.25;

ITM.T =  0.005;
ETM.T =  50E-6;

ITM.L = 50E-6;
ETM.L = 50E-6;

ITM.R = 1 - ITM.T - ITM.L;
ETM.R = 1 - ETM.T - ETM.L;


ITM.r = sqrt(ITM.R);
ETM.r= sqrt(ETM.R);
ITM.t = sqrt(ITM.T);
ETM.t= sqrt(ETM.T);

Refrac_index = 1.45;

% Variable

Laser.k_prop = 2*pi/Laser.lambda;
Grid.step = Grid.Length/Grid.Num_point;


%--------------------- Create the scale for the grid -----------

Grid.vector = 1:Grid.Num_point;

Grid.axis = -Grid.Length/2 + Grid.step/2 + (Grid.vector-1)*Grid.step;
Grid.axis_fft = -1/(2*Grid.step) + (Grid.vector-1)*1/(Grid.Num_point*Grid.step);    

%-------------------- Create mirror mask --------------------------

[Grid.X,Grid.Y] = meshgrid(Grid.axis);
[Grid.FFT_X,Grid.FFT_Y] = meshgrid(Grid.axis_fft);

Grid.D2_square = Grid.X.^2 + Grid.Y.^2;
Grid.D2 = sqrt(Grid.D2_square);

Mirror.mask_index = find (Grid.D2 < (Mirror.Diam/2));
Mirror.mask = zeros(Grid.Num_point,Grid.Num_point,'double');
Mirror.mask(Mirror.mask_index) = 1;

%-------------------- Create propagation matrix  --------------------------
Mat_propagation = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

Mat_propagation = exp(i*(-Laser.k_prop*Length_cav + ...
            pi*Laser.lambda*(Grid.FFT_X.^2 + Grid.FFT_Y.^2)*Length_cav));
Mat_propagation = fftshift(Mat_propagation);
        
%--------------------Create mirrors matrix -------------------------------

CreateMirror;

%------------------- Create input EM field ----------------------

Field.Start = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

Field.Start = Laser.amplitude .* exp(-Grid.D2_square/Laser.beam_radius^2).* ...
            exp(-i*Laser.k_prop.*Grid.D2_square./(2*Laser.wavefront_radius));




