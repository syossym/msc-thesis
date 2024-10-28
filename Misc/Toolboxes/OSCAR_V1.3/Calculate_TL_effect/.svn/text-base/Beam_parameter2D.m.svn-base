function [Beam_radius, Beam_ROfC] = Beam_parameter2D(wave_name)
% Find the beam radius and the radius wavefront of curvature of a Gaussian
% beam.
% !!! This procedure only work for the fundamental TEM_00 Gaussian beam !!!
% New version 2D fit

global Grid;
global Laser;

Amplitude_beam = abs(wave_name);
Phase_beam = angle(wave_name);


Result_fit2 = Fit_Gaussian(Grid.D2_square,Amplitude_beam);

%fprintf('\n Beam waist (m): %g \n',Result_fit2(2));
Beam_radius = Result_fit2(2);
 

% Find on which length do the fit of the wavefront radius of curvature 
% Do the 2D fit between -waist/2 +waist/2

[Val_tmp, Index_tmp_max] = min(abs(Grid.axis - Beam_radius));
[Val_tmp, Index_tmp_min] = min(abs(Grid.axis + Beam_radius));

Phase_cut = Phase_beam(Index_tmp_min:Index_tmp_max,Index_tmp_min:Index_tmp_max);
Phase_r =  Grid.D2(Index_tmp_min:Index_tmp_max,Index_tmp_min:Index_tmp_max);

poly = polyfit(Phase_r,Phase_cut,2);

beam_radius_fit = -Laser.k_prop/(2*poly(1));
% 
%fprintf(' Radius of the beam (m): %g \n',beam_radius_fit);
Beam_ROfC = beam_radius_fit;
end

function Result_fit = Fit_Gaussian(r2_grid,Amplitude_beam)
    global Grid;
    global Laser;
     
    start_point = [Amplitude_beam(Grid.Num_point/2,Grid.Num_point/2);Laser.beam_radius];
    
    options = optimset('TolX',1E-6);
    Result_fit = fminsearch(@gaussfun,start_point,options,r2_grid,Amplitude_beam);

    function sse = gaussfun(params,data_x,data_y)        
        % Function: Amp * exp -(r^2/Width^2)
        % Amp = params(1);
        % Width = params(2);
        
        FittedCurve = params(1) * exp(-data_x/params(2).^2);
        ErrorVector = FittedCurve - data_y;
        sse = sum(sum(ErrorVector .^ 2));
    end
end


