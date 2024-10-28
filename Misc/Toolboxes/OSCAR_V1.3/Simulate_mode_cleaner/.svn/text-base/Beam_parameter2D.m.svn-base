function [Beam_radius, Beam_ROfC] = Beam_parameter(wave_name)
% Find the beam radius and the radius wavefront of curvature of a Gaussian
% beam.
% !!! This procedure only work for the fundamental TEM_00 Gaussian beam !!!

global Grid;
global Laser;


% Find the beam radius first in x and y direction

Power_distri = abs(wave_name).^2;
Power_max = max(max(Power_distri));

Result_fit2 = Fit_Gaussian2D(Power_distri,Power_max);

fprintf('\n Beam waist x (m): %g \n',Result_fit2(3));
% fprintf(' Center of the beam x (m): %g \n',Result_fit2(2));
 fprintf(' Beam waist y (m): %g \n',Result_fit2(5));
%fprintf(' Center of the beam x (m): %g \n',Result_fit2(4));

% find the index center for x and y
temp = abs(Grid.axis - Result_fit2(2));
i_c_x = find(min(temp)==temp, 1 );

temp = abs(Grid.axis - Result_fit2(4));
i_c_y = find(min(temp)==temp, 1 );

% return the average beam radius
Beam_radius = 0.5*(Result_fit2(3)+Result_fit2(5));


% Find the wavefront radius of curvature in x/y

Phase_beam = angle(wave_name);

Cross_sec_phase_x = Phase_beam(i_c_x,:);
Cross_sec_phase_y = Phase_beam(:,i_c_y);

Phase_unwrapped_x = unwrap(Cross_sec_phase_x);
Phase_unwrapped_y = unwrap(Cross_sec_phase_y);

% figure(102)
% plot(Grid.axis,Phase_unwrapped_x-min(Phase_unwrapped_x),Grid.axis,Phase_unwrapped_y-min(Phase_unwrapped_y))

% figure(2)
% plot(Grid.axis,Phase_unwrapped_x)
% 
% figure(3)
% plot(Grid.axis,Phase_unwrapped_y)



% Find on which length do the fit of the wavefront radius of curvature 
% Do the fit between -waist +waist

Index_waist_x = floor(Result_fit2(3)/(Grid.step));
Index_waist_y = floor(Result_fit2(5)/(Grid.step));

for j = 1:(2*Index_waist_x)     
    Phasefit_tmp(j) = Phase_unwrapped_x(i_c_y-Index_waist_x+j);
    Gridfit_tmp(j) = Grid.axis(i_c_y-Index_waist_x+j);
end

poly = polyfit(Gridfit_tmp,Phasefit_tmp,2);
beam_radius_fit_x = -Laser.k_prop/(2*poly(1));

clear Phasefit_tmp
clear Gridfit_tmp

for j = 1:(2*Index_waist_y)    
     Phasefit_tmp(j) = Phase_unwrapped_y(i_c_x-Index_waist_y+j);
     Gridfit_tmp(j) = Grid.axis(i_c_x-Index_waist_y+j);
end

poly = polyfit(Gridfit_tmp,Phasefit_tmp,2);
beam_radius_fit_y = -Laser.k_prop/(2*poly(1));

 fprintf(' Radius of the wavefront x (m): %g \n',beam_radius_fit_x);
 fprintf(' Radius of the wavefront y (m): %g \n',beam_radius_fit_y);

Beam_ROfC = 0.5*(beam_radius_fit_x+beam_radius_fit_y);
end

function Result_fit = Fit_Gaussian2D(Power_distri,Power_max)
    
    global Grid;
    global Laser;
    
    start_point = [Power_max;0;Laser.beam_radius;0;Laser.beam_radius];
    
    options = optimset('TolX',1E-6);
    Result_fit = fminsearch(@gaussfun,start_point,options,Power_distri);

    function sse = gaussfun(params,meas_pow)        
        % Function: B * exp -2(x-x_origin)^2/Width_x^2 * exp-2(y-y_origin)^2/Width_y^2
        % B = params(1);
        % x_origin = params(2);
        % Width_x = params(3);
        % y_origin = param(4);
        % Width_y = param(5);
        
        %global Grid;
        
        FittedCurve = params(1) .* exp(-2*((Grid.X - params(2)).^2)/params(3).^2);
        FittedCurve = FittedCurve .* exp(-2*((Grid.Y - params(4)).^2)/params(5).^2);
        
        ErrorVector = FittedCurve - meas_pow;
        sse = sum(sum(ErrorVector .^ 2));
    end
end


