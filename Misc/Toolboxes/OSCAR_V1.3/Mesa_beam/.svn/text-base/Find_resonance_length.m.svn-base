% Macro to find the resonance length of the cavity
% First we scan the cavity over half a wavelength then we zoom on the maximum
% to find the exact microscopique length to be on resonance

% Phase_shift = exp(sqrt(-1)*k_prop* Dist_reson);

Length.nb_propa_field = 1200;        % Nb of light round trip to take into account
Length.nb_iter = 2000;               % Nb of points used to find the resonance length

fprintf('\n------ Calculating the resonance length ------ \n')

Propagate_Field;

Length.find_reso = zeros(1,Length.nb_iter,'double');
Power.find_reso = zeros(1,Length.nb_iter,'double');
Field.temp = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

% Create the length vector to scan the cavity
Length.find_reso = (1:Length.nb_iter) * Laser.lambda/Length.nb_iter;


fprintf(' Calculating power inside the cavity ...       ')

for q = 1:Length.nb_iter
    Field.temp = Build_Field_Cavity(Length.find_reso(q));
    Power.find_reso(q) = Calculate_power(Field.temp);

    if (rem(q,Length.nb_iter/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*q/(Length.nb_iter))       
    end     
end

fprintf('Finished ... \n')

figure(101)
clf;
semilogy(Length.find_reso,Power.find_reso)
title('Scanning the cavity over one FSR')
xlabel('Microscopique displacement')
ylabel('Circulating power')

tmp_save(:,1) = Length.find_reso';
tmp_save(:,2) = Power.find_reso';
 
save('Cavity_scan.txt','tmp_save','-ASCII');
clear('tmp_save')

% We found the resonance length for the fundamental mesa beam and the
% equivalent LG_10
Detuning.peak1 = 0;
Detuning.peak2 = 2.7664E-7;


    figure(105)
    clf;
    subplot(1,2,1)
        Plot_Field(Detuning.peak1)
        title(['Detuning position: ',num2str(Detuning.peak1)])
    subplot(1,2,2)
        Plot_Field(Detuning.peak2)
        title(['Detuning position: ',num2str(Detuning.peak2)])

  Mode_00 = Build_Field_Cavity(0);       
  Mode_10 = Build_Field_Cavity(2.7664E-7); 
  
  Mode_00 = Mode_00/sqrt(Calculate_power(Mode_00));
  Mode_10 = Mode_10/sqrt(Calculate_power(Mode_10));
  
  figure(106)
  clf;
  plot(Grid.axis,abs(Mode_00(:,Grid.Num_point/2)))
  hold all
  plot(Grid.axis,abs(Mode_10(:,Grid.Num_point/2)))
  hold off
  title('Amplitude profile comparison')
  legend('Fundamental Mesa beam','Mesa Mode LG10')

