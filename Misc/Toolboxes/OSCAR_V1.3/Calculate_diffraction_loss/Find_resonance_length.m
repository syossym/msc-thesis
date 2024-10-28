% Macro to find the resonance length of the cavity
% First we scan the cavity over half a wavelength then we zoom on the maximum
% to find the exact microscopique length to be on resonance

% Phase_shift = exp(sqrt(-1)*k_prop* Dist_reson);

Length.nb_propa_field = 1000;        % Nb of light round trip to take into account
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

% We found several resonant peaks in the cavity circulating power. Let's have a
% look at the 4 highest peaks to find the mode HG10
Detuning.peak1 = 3.181E-007;
Detuning.peak2 = 5.086E-007;
Detuning.peak3 = 6.99E-007;
Detuning.peak4 = 8.89E-007;

    figure(105)
    clf;    
    subplot(2,2,1)
        Plot_Field(Detuning.peak1)
        title(['Detuning position: ',num2str(Detuning.peak1)])
    subplot(2,2,2)
        Plot_Field(Detuning.peak2)
        title(['Detuning position: ',num2str(Detuning.peak2)])
    subplot(2,2,3)
        Plot_Field(Detuning.peak3)
        title(['Detuning position: ',num2str(Detuning.peak3)])
    subplot(2,2,4)
        Plot_Field(Detuning.peak4)
        title(['Detuning position: ',num2str(Detuning.peak4)])

% So we manually found the mode HG10 is for a detuning of 3.181E-007 
        
%----------------------------Zoom on the peak --------------------------

% Zoom on the LG10 resonance peak
Length.resonance =  3.181E-007; 


Length.Span_zoom = 1e-09;
Length.nb_iter_zoom = 1000;
% 
Length.find_reso_zoom = zeros(1,Length.nb_iter_zoom,'double');
Power.find_reso_zoom = zeros(1,Length.nb_iter_zoom,'double');


% Create the length vector to scan the cavity
Length.find_reso_zoom = Length.resonance - Length.Span_zoom/2 + (1:Length.nb_iter_zoom)*Length.Span_zoom/Length.nb_iter_zoom;

fprintf(' Zooming on the resonance peak ...       ')

for q = 1:Length.nb_iter_zoom
    Field.temp = Build_Field_Cavity(Length.find_reso_zoom(q));
    Power.find_reso_zoom(q) = Calculate_power(Field.temp);

    if (rem(q,Length.nb_iter_zoom/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*q/(Length.nb_iter_zoom))       
    end     
end

fprintf('Finished ... \n')
figure(102)
plot(Length.find_reso_zoom,Power.find_reso_zoom)
title('zooming on the resonance peak')
xlabel('Microscopique displacement')
ylabel('Circulating power')


[Power.max_reso_zoom,Length.index_reso_zoom] = max(Power.find_reso_zoom);
Length.reso_zoom = Length.find_reso_zoom(Length.index_reso_zoom);

fprintf(' Microscopique resonance length: %-9.7d ... \n',Length.reso_zoom);

