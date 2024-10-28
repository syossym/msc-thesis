function Result = Make_propagation(Mat_name,Distance)

global Grid
global Laser

Mat_propagation = exp(i*(-Laser.k_prop*Distance + ...
            pi*Laser.lambda*(Grid.FFT_X.^2 + Grid.FFT_Y.^2)*Distance));

Wave_fft = fftshift(fft2(Mat_name));
Wave_prop = Mat_propagation .* Wave_fft;
Result = ifft2(ifftshift(Wave_prop));
