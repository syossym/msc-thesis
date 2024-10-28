% My_FFT_code.m
% Simple code to understand how to progate a laser beam over a certain distance 
% using a FFT method
% !! For clarity, this code is absolutely not optmised !!
% Jerome Degallaix - March 2008

clear all

%-----------------------------------------------------------------
fprintf('\n*******************************************************\n ')
fprintf('                   Simple FFT code                       \n ')

%--------------------- Parameters simulations --------------------
% Dimension of the 2D grid

Grid.Num_point = 512;                       % Number of point of the grid
Grid.Length = 0.16;                         % Physical dimension of the grid in meter

% Parameters of the input laser beam
Laser.lambda = 1064e-9;                     % Wavelength of the laser in meter

% Distance to propagate the beam
Distance_prop = 100;                        % The propagation will be simulated over this distance (m)

% Additional variables
Laser.k_prop = 2*pi/Laser.lambda;           % Propagation constant
Grid.step = Grid.Length/Grid.Num_point;     % Physical size of one pixel of the grid


%--------------------- Create the scale for the grid -----------

Grid.vector = 1:Grid.Num_point;             

Grid.axis = -Grid.Length/2 + Grid.step/2 + (Grid.vector-1)*Grid.step;               % Physical position of the pixels along the x/y axis
Grid.axis_fft = -1/(2*Grid.step) + (Grid.vector-1)*1/(Grid.Num_point*Grid.step);    % Spatial frequency of the pixels in the Fourrier plane

%-------------------- Create propagation matrix  --------------------------

Mat_propagation = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

for m = 1:Grid.Num_point
    for n = 1:Grid.Num_point

        Mat_propagation(m,n) = exp(i*(-Laser.k_prop*Distance_prop + ...
            pi*Laser.lambda*(Grid.axis_fft(m)^2 + Grid.axis_fft(n).^2)*Distance_prop));

    end
end

%------------------- Create the input field ----------------------
% The input field is a uniform square of light of 10cm * 10cm
Square_dimension = 0.04;              % Physical dimension of the side of the square light field 

Field.Start = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Laser.amplitude = 1;        % Arbitrary amplitude

for m = 1:Grid.Num_point
    for n = 1:Grid.Num_point

        if (abs(Grid.axis(m)) < Square_dimension/2) & (abs(Grid.axis(n)) < Square_dimension/2)
            Field.Start(m,n) = 1;
        end
                    
    end
end

% % Try also for a round uniform light intensity
% 
% Round_diameter = 0.1;              
% Field.Start = zeros(Grid.Num_point,Grid.Num_point,'double');
% Laser.amplitude = 1;        
% 
% for m = 1:Grid.Num_point
%     for n = 1:Grid.Num_point
%         if ((Grid.axis(m)^2 + Grid.axis(n)^2 ) < (Round_diameter/2)^2)
%             Field.Start(m,n) = 1;
%         end                    
%     end
% end

%------------------ Propagate the field ---------------------------

Field.Fourier = fftshift(fft2(Field.Start));       % Do the Fourier transform of the input field
Field.Fourier = Field.Fourier .* Mat_propagation;  % Do the propagation in the frequency domain
Field.End = ifft2(ifftshift(Field.Fourier));       % Do the inverse Fourier transform 


%----------------- Display the results -------------------------------

figure(101)
clf;
imagesc(Grid.axis,Grid.axis,Field.Start)
title('Input field')
axis square

figure(102)
clf;
imagesc(Grid.axis,Grid.axis,abs(Field.End))
title('Propagated field')
axis square

