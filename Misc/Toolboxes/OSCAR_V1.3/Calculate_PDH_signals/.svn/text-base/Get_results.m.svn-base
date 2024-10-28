Iter.final = 400;               % Number of round trip for each detuning
Iter.detuning_num = 200;        % Number of length to scan
% Scan the micoscopic length of the cavity over one wavelength
% to have a better plot, the scan is centered on the carrier resonance frequency:

Iter.detuning_length = (8.778E-7-Laser.lambda/2)+(1:Iter.detuning_num) * Laser.lambda/Iter.detuning_num;

Field.Circ = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Total = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Reflect = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Leak = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.Transmit = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

Power.circ = zeros(1,Iter.detuning_num,'double');

% For the sidebands
Field.Total_sidebands_plus = zeros(Grid.Num_point,Grid.Num_point,'double');
Field.Total_sidebands_minus = zeros(Grid.Num_point,Grid.Num_point,'double');

Field.Reflect_sidebands_plus = zeros(Grid.Num_point,Grid.Num_point,'double');
Field.Reflect_sidebands_minus = zeros(Grid.Num_point,Grid.Num_point,'double');

% Additional phase shift induced by the sidebands during one round trip
Phase_shift_sidebands_plus = exp(i*Sidebands.k_prop* 2*Length_cav);
Phase_shift_sidebands_minus = exp(-i*Sidebands.k_prop* 2*Length_cav);


fprintf('\n ------ Scanning the cavity  ------ \n  ')
fprintf('Be patient ... calculating ...       ')

for r = 1:Iter.detuning_num
    % Phase shift per round trip for the carrier
    Phase_shift =  exp(i*Laser.k_prop* Iter.detuning_length(r));

    Field.Total = zeros(Grid.Num_point,Grid.Num_point,'double');
    Field.Reflect = zeros(Grid.Num_point,Grid.Num_point,'double');
    Field.Leak = zeros(Grid.Num_point,Grid.Num_point,'double');

    % For the sidebands
    Field.Total_sidebands_plus = zeros(Grid.Num_point,Grid.Num_point,'double');
    Field.Total_sidebands_minus = zeros(Grid.Num_point,Grid.Num_point,'double');

    Field.Reflect_sidebands_plus = zeros(Grid.Num_point,Grid.Num_point,'double');
    Field.Reflect_sidebands_minus = zeros(Grid.Num_point,Grid.Num_point,'double');


    Field.Circ = Propa_mirror(Field.Start, Mirror.ITM_trans,i*ITM.t);
    Field.Reflect = Propa_mirror(Field.Start,Mirror.ITM_ref,ITM.r);

    Field.Circ_sidebands_plus = Field.Circ*Ration_amp;   
    Field.Circ_sidebands_minus = -Field.Circ*Ration_amp;
    
    for q = 1:Iter.final

        Field.Total = Field.Total + Field.Circ;
        Field.Total_sidebands_plus = Field.Total_sidebands_plus + Field.Circ_sidebands_plus;
        Field.Total_sidebands_minus = Field.Total_sidebands_minus + Field.Circ_sidebands_minus;

        Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
        Field.Circ = Propa_mirror(Field.Circ,Mirror.ETM_cav,ETM.r);
        Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
        
        Field.Circ = Propa_mirror(Field.Circ, Mirror.ITM_cav,ITM.r);

        % Add the proper round trip phase shift
        Field.Circ = Field.Circ * Phase_shift;
        Field.Circ_sidebands_plus = Field.Circ*Ration_amp*(Phase_shift_sidebands_plus).^q;
        Field.Circ_sidebands_minus = -Field.Circ*Ration_amp*(Phase_shift_sidebands_minus).^q;
        
        
    end



    %-------------------------------------------------------------------
    % Calculate the transmitted and reflected field for the carrier

    Field.Leak = Make_propagation(Field.Total,Mat_propagation);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ETM_cav,ETM.r);
    Field.Leak = Make_propagation(Field.Leak,Mat_propagation);
    Field.Leak = Field.Leak * Phase_shift;
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ITM_trans,i*ITM.t);
    Field.Reflect_carrier = Field.Reflect + Field.Leak;


    % Do the same for the upper and lower sidebands

    Field.Leak = Make_propagation(Field.Total_sidebands_plus,Mat_propagation);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ETM_cav,ETM.r);
    Field.Leak = Make_propagation(Field.Leak,Mat_propagation);
    Field.Leak = Field.Leak * Phase_shift* Phase_shift_sidebands_plus;
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ITM_trans,i*ITM.t);
    Field.Reflect_sidebands_plus = Ration_amp*Field.Reflect + Field.Leak;

    Field.Leak = Make_propagation(Field.Total_sidebands_minus,Mat_propagation);
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ETM_cav,ETM.r);
    Field.Leak = Make_propagation(Field.Leak,Mat_propagation);
    Field.Leak = Field.Leak * Phase_shift* Phase_shift_sidebands_minus;
    Field.Leak = Propa_mirror(Field.Leak,Mirror.ITM_trans,i*ITM.t);
    Field.Reflect_sidebands_minus = -Ration_amp*Field.Reflect + Field.Leak;

    
    % Calculate PDH signal for this detuning

    Power.total(r) = Calculate_power(Field.Total)+ Calculate_power(Field.Total_sidebands_plus)+Calculate_power(Field.Total_sidebands_minus);

    Field.PDH = imag(Field.Reflect_carrier.*conj(Field.Reflect_sidebands_plus)+conj(Field.Reflect_carrier).*Field.Reflect_sidebands_minus);
    Power.PDH(r) = sum(sum(Field.PDH))/(Laser.amplitude*Sidebands.amplitude);
     
   if (rem(r,Iter.detuning_num/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b  %-3.0i %% ',100*r/Iter.detuning_num)
    end
    
end

fprintf('\b\b\b\b\b\b\b\b  Finished! ')
fprintf('\n')
    
     
%-------------------------------------------------------------------
 
figure(105)
plot(Iter.detuning_length,Power.total)
title('Circulating power inside the cavity')
xlabel('Detuning (nm)')
ylabel('Circulating power (a.u.)')

figure(106)
plot(Iter.detuning_length,Power.PDH)
title('PDH error signal as a function of the detuning')
xlabel('Detuning (nm)')
ylabel('PDH signal (a.u.)')