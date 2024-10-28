%Create mirror profile for flat beam

Field.mesa_tmp = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ITM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Mirror.ETM_cav = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));


%---------------------- Create nearly concentric flat beam -------------------
% Create the profile at the cavity waist

waist_0 = sqrt((Laser.lambda  * Length_cav) / (2*pi));
p = 3*waist_0; % Defined the width of the Mesa beam

x_temp = Grid.D2/waist_0;
Field.mesa_tmp = (1./x_temp).*exp(-x_temp.^2).*besselj(1,2*x_temp*p/waist_0);

clear('x_temp')

% Propagate the beam from the waist to the mirror
Field.On_mirror = Make_propagation(Field.mesa_tmp,Mat_propagation_h);

% Take the phase of the beam at the mirror and calculate the equivalent
% change in sagitta
Mirror.ETM_cav = (2/Laser.k_prop)*angle(Field.On_mirror);
Mirror.ITM_cav = Mirror.ETM_cav;

clear('Field.On_mirror','Field.mesa_tmp');




