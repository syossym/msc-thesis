function Result = Make_propagation(Mat_name,Mat_propagation)

% Wave_fft = fftshift(fft2(Mat_name));
% Wave_prop = Mat_propagation .* Wave_fft;
% Result = ifft2(ifftshift(Wave_prop));


Wave_fft = fft2(Mat_name);
Wave_prop = Mat_propagation .* Wave_fft;
Result = ifft2(Wave_prop);