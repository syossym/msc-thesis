Iter.final = 1000;

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

Field.Circ = Propa_mirror(Field.Start, Mirror.M1_trans,i*M1.t,Grid_beta.X);
Field.Reflect = Propa_mirror(Field.Start,Mirror.M1_ref,M1.r,Grid_beta.X);

for q = 1:Iter.final

     Field.Total = Field.Total + Field.Circ;
     Power.Buildup(q) = Calculate_power(Field.Total);
     
     Field.Circ = Make_propagation(Field.Circ,MC.Length_small);
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M2_cav,M2.r,Grid_beta.X);
     Field.Circ = Make_propagation(Field.Circ,MC.Length_side);
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M3_cav,M3.r,Grid_alpha.X);
     Field.Circ = Make_propagation(Field.Circ,MC.Length_side);
     
     Field.Circ = Field.Circ * Phase_shift;
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M1_cav,M1.r,Grid_beta.X);
    
    if (rem(q,Iter.final/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b  %-3.0i %% ',100*q/Iter.final)
    end

end
    
fprintf('\b\b\b\b\b\b\b\b  Finished! ')  
fprintf('\n')

%-------------------------------------------------------------------
    % Calculate the transmitted and reflected field

    Field.Transmit = Make_propagation(Field.Total,MC.Length_small);
    Field.Transmit = Propa_mirror(Field.Transmit,Mirror.M2_trans,i*M2.t,Grid_beta.X);

    Field.Leak = Make_propagation(Field.Total,MC.Length_small);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.M2_cav,M2.r,Grid_beta.X);
    Field.Leak = Make_propagation(Field.Leak,MC.Length_side);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.M3_cav,M3.r,Grid_alpha.X);
    Field.Leak = Make_propagation(Field.Leak,MC.Length_side);
     
    Field.Leak = Field.Leak * Phase_shift;
    Field.Leak = Propa_mirror(Field.Leak,Mirror.M1_cav,i*M1.t,Grid_beta.X);  
    
    Field.Reflect = Field.Reflect + Field.Leak;
    
%-------------------------------------------------------------------
    % Write the total circulating power
    % Write the beam radius on the ITM
    % and also the cavity waist size and position

    fprintf('\n ------ Display the results ------ \n')
% 
    fprintf(' Circulating power (W): %f \n',Calculate_power(Field.Total));
    fprintf(' Reflected power (W): %f \n',Calculate_power(Field.Reflect)); 
    fprintf(' Transmitted power (W): %f \n',Calculate_power(Field.Transmit));
%     

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

      
        

