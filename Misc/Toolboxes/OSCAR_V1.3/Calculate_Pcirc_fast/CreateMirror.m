%Create mirror grid

Mirror.ITM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ETM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ITM_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ETM_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ITM_ref = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

% Optical path difference of the mirrors in reflection, inside the cavity

Mirror.ITM_cav = -(ITM.RofC - sign(ITM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square))*2;
Mirror.ETM_cav = -(ETM.RofC - sign(ETM.RofC)*sqrt(ETM.RofC^2 - Grid.D2_square))*2;


% Optical path difference of the mirrors in in transmission

Mirror.ETM_trans = (Refrac_index-1)*(ETM.RofC - sign(ITM.RofC)*sqrt(ETM.RofC^2 - Grid.D2_square));
Mirror.ITM_trans = (Refrac_index-1)*(ITM.RofC - sign(ETM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square));


% Optical path difference of the input mirror in reflection in the
% substrate

Mirror.ITM_ref = (Refrac_index)*(ITM.RofC - sign(ITM.RofC)*sqrt(ITM.RofC^2 - Grid.D2_square))*2;

% To make the calculation faster, calculate the wavefront phase shit due to
% reflection/transmission

Mirror.ITM_cav = exp(-i * Mirror.ITM_cav*Laser.k_prop).* Mirror.mask;
Mirror.ETM_cav = exp(-i * Mirror.ETM_cav*Laser.k_prop).* Mirror.mask;

Mirror.ITM_trans = exp(-i * Mirror.ITM_trans*Laser.k_prop).* Mirror.mask;
Mirror.ETM_trans = exp(-i * Mirror.ETM_trans*Laser.k_prop).* Mirror.mask;

Mirror.ITM_ref = exp(-i * Mirror.ITM_ref*Laser.k_prop).* Mirror.mask;

