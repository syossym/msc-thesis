%Create mirror grid

Mirror.ITM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ETM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ITM_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ETM_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ITM_ref = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

% Optical path difference of the cold mirrors in reflection, inside the cavity

Mirror.ITM_cav = -(ITM.RofC - sign(ITM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square))*2;
Mirror.ETM_cav = -(ETM.RofC - sign(ETM.RofC)*sqrt(ETM.RofC^2 - Grid.D2_square))*2;


% Optical path difference of the mirrors in in transmission

Mirror.ETM_trans = (Refrac_index-1)*(ETM.RofC - sign(ITM.RofC)*sqrt(ETM.RofC^2 - Grid.D2_square));
Mirror.ITM_trans = (Refrac_index-1)*(ITM.RofC - sign(ETM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square));


% Optical path difference of the input mirror in reflection in the
% substrate

Mirror.ITM_ref = (Refrac_index)*(ITM.RofC - sign(ITM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square))*2;

% Now we include thermal lensing. We load the file with the optical path
% difference of the input mirror in transmission as well as the sagitta
% change of the mirrors (see manual for details)

load('From_ANSYS.txt')
loaded.radius = From_ANSYS(:,1);
loaded.TL = interp1(From_ANSYS(:,1),From_ANSYS(:,2),Grid.D2,'spline')*1;
loaded.sag = interp1(From_ANSYS(:,1),From_ANSYS(:,3),Grid.D2,'spline')*1;

% Add the thermal lensing distortion to the previous wavefront
% Special attention to the sign!

Mirror.ITM_cav = Mirror.ITM_cav - 2*loaded.sag;
Mirror.ETM_cav = Mirror.ETM_cav - 2*loaded.sag;

Mirror.ETM_trans = Mirror.ETM_trans + loaded.TL - (Refrac_index-1)*loaded.sag;
Mirror.ITM_trans = Mirror.ITM_trans + loaded.TL - (Refrac_index-1)*loaded.sag;

Mirror.ITM_ref = Mirror.ITM_ref + 2*loaded.TL - 2*Refrac_index*loaded.sag;


