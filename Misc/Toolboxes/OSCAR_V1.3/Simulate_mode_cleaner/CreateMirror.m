%Create mirror grid

Mirror.M1_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M2_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M3_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M1_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M2_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M3_trans = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.M1_ref = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

% Optical path difference of the mirrors in reflection, inside the cavity

Mirror.M1_cav = -(M1.RofC - sign(M1.RofC)*sqrt(M1.RofC^2 - Grid.D2_square))*2;
Mirror.M2_cav = -(M2.RofC - sign(M2.RofC)*sqrt(M2.RofC^2 - Grid.D2_square))*2;
Mirror.M3_cav = -(M3.RofC - sign(M3.RofC)*sqrt(M3.RofC^2 - Grid.D2_square))*2;

% Optical path difference of the mirrors in in transmission

Mirror.M1_trans = (Refrac_index-1)*(M1.RofC - sign(M1.RofC)*sqrt(M1.RofC^2 - Grid.D2_square));
Mirror.M2_trans = (Refrac_index-1)*(M2.RofC - sign(M2.RofC)*sqrt(M2.RofC^2 - Grid.D2_square));
Mirror.M3_trans = (Refrac_index-1)*(M3.RofC - sign(M3.RofC)*sqrt(M3.RofC^2 - Grid.D2_square));

% Optical path difference of the input mirror in reflection in the
% substrate

Mirror.M1_ref = (Refrac_index)*(M1.RofC - sign(M1.RofC)*sqrt(M1.RofC^2 - Grid.D2_square))*2;


