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

Laser.lambda = 1064e-9;

% The parameters for the input laser beam are defined at the end of the
% file
%--------------------- Cavity parameters ---------------------------

Length_cav = 1000;

ITM.RofC = 510;
ETM.RofC = 510;

Mirror.Diam = 0.25;

ITM.T =  0.0005;
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
        
%--------------------Create mirrors matrix -------------------------------

CreateMirror;

%------------------- Create input EM field ----------------------

Field.Start = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

% Input field used to excite the HG10 mode

thick_light = 0.15;
amplitude_E = 1;

tmp_index = find (abs(Grid.X) < (thick_light/2) );
Field.Start(tmp_index) = sign(Grid.X(tmp_index)).*amplitude_E;
clear('tmp_index')

% Also to try a uniform circular field to excite the LGm0 modes

% diam_light = 0.10;
% amplitude_E = 1;
% 
% tmp_index = find (Grid.D2 < (diam_light/2) );
% Field.Start(tmp_index) = amplitude_E;
% clear('tmp_index')













%     case '1'
%             % First distribution: circular symmetry (to excite LG modes)
%             diam_light = 0.20;
%             amplitude_E = 1;
% 
%             for j = 1:Num_point
%                 for k = 1:Num_point
%                     Rad(j,k) = sqrt(Grid_axis(j)^2 + Grid_axis(k)^2);
%                     if Rad(j,k) < (diam_light/2)            
%                         Field_Start(j,k) = amplitude_E;
%                     else
%                         Field_Start(j,k) = 0;
%                     end
%                 end
%             end
%     case '2'
%             % Second distribution: vertical strip of light (to excite HG10
%             % modes)





