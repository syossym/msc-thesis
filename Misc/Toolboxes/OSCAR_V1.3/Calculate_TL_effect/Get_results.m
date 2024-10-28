
CreateField;
Find_resonance_length;

%Length.reso_zoom =  9.0985000e-007;
%Length.reso_zoom = 8.9373600e-007;
Iter.final = 5000;

Field.Circ = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Total = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Reflect = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Leak = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Transmit = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

Power.Buildup = zeros(1,Iter.final,'double');


% Phase shift per round trip to make it resonant in the cavity
Phase_shift =  exp(i*Laser.k_prop* Length.reso_zoom);

fprintf('\n ------ Calculating circulating field inside the cavity ------ \n  ')
fprintf('Be patient ... calculating ...       ')

Field.Circ = Propa_mirror(Field.Start, Mirror.ITM_trans,i*ITM.t);
Field.Reflect = Propa_mirror(Field.Start,Mirror.ITM_ref,ITM.r);

for q = 1:Iter.final

    Field.Total = Field.Total + Field.Circ;
    Power.Buildup(q) = Calculate_power(Field.Total);
    
    Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
    Field.Circ = Propa_mirror(Field.Circ,Mirror.ETM_cav,ETM.r);
    Field.Circ = Make_propagation(Field.Circ,Mat_propagation);

    Field.Circ = Field.Circ * Phase_shift;   
    %Field.Reflect = Field.Reflect + Propa_mirror(Field.Circ,Mirror.ITM_trans,i*ITM.t);
    
    Field.Circ = Propa_mirror(Field.Circ, Mirror.ITM_cav,ITM.r);   
    
    if (rem(q,Iter.final/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b  %-3.0i %% ',100*q/Iter.final)
    end

end
    
fprintf('\b\b\b\b\b\b\b\b  Finished! ')  
fprintf('\n')

%-------------------------------------------------------------------
    % Calculate the transmitted and reflected field

    Field.Transmit = Make_propagation(Field.Total,Mat_propagation);
    Field.Transmit = Propa_mirror(Field.Transmit,Mirror.ETM_trans,i*ETM.t);

    Field.Leak = Make_propagation(Field.Total,Mat_propagation);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ETM_cav,ETM.r);
    Field.Leak = Make_propagation(Field.Leak,Mat_propagation);
    Field.Leak = Field.Leak * Phase_shift;
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ITM_trans,i*ITM.t);
    Field.Reflect = Field.Reflect + Field.Leak;
    
%-------------------------------------------------------------------
    % Write the total circulating power
    % Write the beam radius on the ITM
    % and also the cavity waist size and position

    fprintf('\n ------ Display the results ------ \n')

    fprintf(' Circulating power (W): %f \n',Calculate_power(Field.Total));
    fprintf(' Reflected power (W): %f \n',Calculate_power(Field.Reflect)); 
    fprintf(' Transmitted power (W): %f \n',Calculate_power(Field.Transmit));
    
    [Field_cav.w, Field_cav.R] = Beam_parameter(Field.Total);
    fprintf(' Beam radius on the ITM (m): %f \n',Field_cav.w);
    
    [Field_cav.w2, Field_cav.R2] = Beam_parameter(Make_propagation(Field.Total,Mat_propagation));
    fprintf(' Beam radius on the ETM (m): %f \n',Field_cav.w2);
    
    Field_cav.q = 1/(1/Field_cav.R - i*Laser.lambda/(pi*Field_cav.w^2));

    fprintf(' Cavity waist size (m): %f \n',sqrt(( (imag(Field_cav.q))*Laser.lambda/pi )));
    fprintf(' Location of the waist from ITM (m): %f \n',real(Field_cav.q));
    
    figure(104)
    clf;    
    subplot(2,2,1)
        Plot_Field(Field.Start)
        title('Input field')
    subplot(2,2,2)
        Plot_Field(Field.Total)
        title('Circulating field')
    subplot(2,2,3)
        Plot_Field(Field.Reflect)
        title('Reflected field')
    subplot(2,2,4)
        Plot_Field(Field.Transmit)
        title('Transmitted field')
  
        
% Field.loss = Field.Total;
% Field.loss = Field.loss/sqrt(Calculate_power(Field.loss));
% 
% Field.loss = Make_propagation(Field.loss,Mat_propagation);
% Field.loss = Propa_mirror(Field.loss,Mirror.ITM_cav,1);
% Field.loss = Make_propagation(Field.loss,Mat_propagation);
% Field.loss = Propa_mirror(Field.loss, Mirror.ITM_cav,1);   
%     
% Dif_loss = (1 - Calculate_power(Field.loss));
% fprintf('Diffraction loss per round trip: %d \n',Dif_loss);    
%     
% 
% 
% [Field_cav.w, Field_cav.R] = Beam_parameter(Field.Start);
% %     fprintf(' Beam radius on the ITM (m): %f \n',Field_cav.w);
% %     
%     Field_cav.q = 1/(1/Field_cav.R - i*Laser.lambda/(pi*Field_cav.w^2));
% 
%     w0 = sqrt(( (imag(Field_cav.q))*Laser.lambda/pi ))
%     Zr = imag(Field_cav.q)
%     Z = real(Field_cav.q)
%     
%     fprintf(' Cavity waist size (m): %f \n',sqrt(( (imag(Field_cav.q))*Laser.lambda/pi )));
%     fprintf(' Location of the waist from ITM (m): %f \n',real(Field_cav.q));

