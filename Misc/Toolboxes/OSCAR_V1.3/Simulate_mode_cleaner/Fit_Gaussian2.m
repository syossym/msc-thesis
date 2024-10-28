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
