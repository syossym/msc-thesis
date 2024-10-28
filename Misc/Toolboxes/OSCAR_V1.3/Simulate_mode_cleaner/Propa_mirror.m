function Output = Propa_mirror(Wave_field, Wave_mirror, reflec,Grid_angle_X)

global Mirror
global Laser
global Grid

% Stretch the laser beam as seen by the mirror
Output = interp2(Grid_angle_X,Grid.Y,Wave_field,Grid.X,Grid.Y,'cubic');

Output = Output .* exp(-i * Wave_mirror*Laser.k_prop) * reflec .* Mirror.mask;


% Go back to the normal grid
Output = interp2(Grid.X,Grid.Y,Output,Grid_angle_X,Grid.Y,'cubic',0);

% Flip the matrix along the x axis
Output = Output(:,Grid.Num_point:-1:1);



