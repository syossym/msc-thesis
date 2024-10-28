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
    % Calculate the diffraction loss
        
Field.loss = Field.Total;
Field.loss = Field.loss/sqrt(Calculate_power(Field.loss));

Field.loss = Make_propagation(Field.loss,Mat_propagation);
Field.loss = Propa_mirror(Field.loss,Mirror.ETM_cav,1);
Field.loss = Make_propagation(Field.loss,Mat_propagation);
Field.loss = Propa_mirror(Field.loss, Mirror.ITM_cav,1);   
    
Dif_loss = (1 - Calculate_power(Field.loss));
fprintf('Diffraction loss per round trip: %d \n',Dif_loss);    

figure(110)
clf;
Plot_Field(Field.loss)

%     
