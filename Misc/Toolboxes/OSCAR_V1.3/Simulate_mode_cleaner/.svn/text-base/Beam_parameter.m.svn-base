function [Beam_radius, Beam_ROfC] = Beam_parameter(wave_name)
% Find the beam radius and the radius wavefront of curvature of a Gaussian
% beam.
% !!! This procedure only work for the fundamental TEM_00 Gaussian beam !!!

global Grid;
global Laser;

Power_temp = Calculate_power(wave_name);
Amplitude_beam = abs(wave_name);
Phase_beam = angle(wave_name);

Cross_sec_Amplitude = Amplitude_beam(Grid.Num_point/2,:);
Cross_sec_phase = Phase_beam(Grid.Num_point/2,:);

Result_fit2 = Fit_Gaussian(Grid.axis,Cross_sec_Amplitude,Power_temp);

%fprintf('\n Beam waist (m): %g \n',Result_fit2(3));
Beam_radius = Result_fit2(3);
 
Phase_unwrapped = unwrap(Cross_sec_phase);

% Find on which length do the fit of the wavefront radius of curvature 
% Do the fit between -waist/2 +waist/2

[Val_tmp, Index_tmp_max] = min(abs(Grid.axis - Result_fit2(3)));
[Val_tmp, Index_tmp_min] = min(abs(Grid.axis + Result_fit2(3)));

Length_interval = (Index_tmp_max - Index_tmp_min);

for j = 1:Length_interval    
     Phasefit_tmp(j) = Phase_unwrapped(Index_tmp_min+j);
     Gridfit_tmp(j) = Grid.axis(Index_tmp_min+j);
end
 
poly = polyfit(Gridfit_tmp,Phasefit_tmp,2);
beam_radius_fit = -Laser.k_prop/(2*poly(1));

%fprintf(' Radius of the beam (m): %g \n',beam_radius_fit);
Beam_ROfC = beam_radius_fit;
end

function Result_fit = Fit_Gaussian(x_axis,Cross_power,Power_temp)
    
    global Laser;
     
    start_point = [sqrt(Power_temp/(pi*Laser.beam_radius^2));0;Laser.beam_radius];
    
    options = optimset('TolX',1E-6);
    Result_fit = fminsearch(@gaussfun,start_point,options,x_axis,Cross_power);

    function sse = gaussfun(params,data_x,data_y)        
        % Function: B * exp -(x-x_origin)^2/Width^2
        % B = params(1);
        % x_origin = params(2);
        % Width = params(3);
        
        FittedCurve = params(1) .* exp(-((data_x - params(2)).^2)/params(3).^2);
        ErrorVector = FittedCurve - data_y;
        sse = sum(ErrorVector .^ 2);
    end
end


