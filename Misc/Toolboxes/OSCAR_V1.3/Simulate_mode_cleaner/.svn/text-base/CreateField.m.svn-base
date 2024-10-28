clear all

global Grid
global Laser
global Mirror
global Length;
global Field
global Grid_alpha 
global Grid_beta

%-----------------------------------------------------------------
fprintf('\n*******************************************************\n ')
fprintf('                   OSCAR FFT code                       \n ')

%--------------------- Parameters simulations --------------------

Grid.Num_point = 128;
Grid.Length = 0.004;

% Parameters of the input laser beam
Laser.lambda = 1064e-9;
Laser.beam_radius = 0.0004;
Laser.wavefront_radius = -5;

%--------------------- 3 mirrors ring cavity ---------------------------
% 2 flat mirrors, 1 curved mirror
% Length_small = distance between the 2 flat mirrors (M1 and M2)
% Length_long = distance between the curved mirror (M3) and the segment between
% the 2 flat mirrors

MC.Length_small = 0.06;
MC.Length_long = 0.20;

M3.RofC = 1;
M1.RofC = 1E99;
M2.RofC = 1E99;

Mirror.Diam = 0.10;

M1.T =  0.01;
M2.T =  0.01;
M3.T =  0.00;

M1.L =  0.00;
M2.L =  0.00;
M3.L =  0.00;

M1.R = 1 - M1.T - M1.L;
M2.R = 1 - M2.T - M2.L;
M3.R = 1 - M3.T - M3.L;

M1.r = sqrt(M1.R);
M2.r = sqrt(M2.R);
M3.r = sqrt(M3.R);

M1.t = sqrt(M1.T);
M2.t = sqrt(M2.T);
M3.t = sqrt(M3.T);

Refrac_index = 1.45;

% Variable

Laser.k_prop = 2*pi/Laser.lambda;
Grid.step = Grid.Length/Grid.Num_point;

MC.alpha = atan(MC.Length_small/(2*MC.Length_long));
MC.beta = 0.5*(pi/2 - MC.alpha);
MC.Length_side = sqrt((MC.Length_small/2)^2 + MC.Length_long^2);

%--------------------- Create the scale for the grid -----------

Grid.vector = 1:Grid.Num_point;

Grid.axis = -Grid.Length/2 + Grid.step/2 + (Grid.vector-1)*Grid.step;
Grid.axis_fft = -1/(2*Grid.step) + (Grid.vector-1)*1/(Grid.Num_point*Grid.step);    

%-------------------- Create mirror mask --------------------------

[Grid.X,Grid.Y] = meshgrid(Grid.axis);
[Grid.FFT_X,Grid.FFT_Y] = meshgrid(Grid.axis_fft);

Grid.D2_square = Grid.X.^2 + Grid.Y.^2;
Grid.D2 = sqrt(Grid.D2_square);
Grid.angle = atan2(Grid.Y,Grid.X);


Mirror.mask_index = find (Grid.D2 < (Mirror.Diam/2));
Mirror.mask = zeros(Grid.Num_point,Grid.Num_point,'double');
Mirror.mask(Mirror.mask_index) = 1;

%------------------- Create the 'stretch' matrix ------------------
% Used to simulate how the beam looks like on the mirror when coming with
% an angle
Grid_alpha.X = Grid.X/cos(MC.alpha);
Grid_alpha.D2_square = Grid_alpha.X.^2 + Grid.Y.^2;
Grid_alpha.D2 = sqrt(Grid_alpha.D2_square);

Grid_beta.X = Grid.X/cos(MC.beta);
Grid_beta.D2_square = Grid_beta.X.^2 + Grid.Y.^2;
Grid_beta.D2 = sqrt(Grid_beta.D2_square);


%--------------------Create mirrors matrix -------------------------------

CreateMirror;

%------------------- Create input EM field ----------------------

Field.Start = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

%Laser.power = 1; 
%Laser.amplitude = sqrt(Laser.power*2/(pi*Laser.beam_radius^2)); % Amplitude for TEM00 input beam
%Field.Start = Laser.amplitude .* exp(-Grid.D2_square/Laser.beam_radius^2).* ...
%            exp(-i*Laser.k_prop.*Grid.D2_square./(2*Laser.wavefront_radius));
        
% Create HG10        
Field.Start = (sqrt(2)*Grid.X/Laser.beam_radius)... 
            .* exp(-Grid.D2_square/Laser.beam_radius^2).* ...
            exp(-i*Laser.k_prop.*Grid.D2_square./(2*Laser.wavefront_radius));

Field.Start = Field.Start + (sqrt(2)*Grid.Y/Laser.beam_radius)... 
            .* exp(-Grid.D2_square/Laser.beam_radius^2).* ...
            exp(-i*Laser.k_prop.*Grid.D2_square./(2*Laser.wavefront_radius));

% Normalise to 1W        
Field.Start = Field.Start/sqrt(Calculate_power(Field.Start));


