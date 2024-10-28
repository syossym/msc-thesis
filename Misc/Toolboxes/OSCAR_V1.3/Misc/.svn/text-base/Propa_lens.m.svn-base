function Output = Propa_lens(Wave_field,f_lens)
% Prograte the input beam 'Wave_field' through a thin lens of focal length
% 'f_lens'. The lens diameter is defined by default as the mirror diameter

global Mirror
global Laser
global Grid

% Defined the wavefront induced by the lens
Wave_lens = -((2*f_lens) - sign(f_lens)*sqrt((2*f_lens)^2 - Grid.D2_square))*2;

Output = Wave_field .* exp(-1i * Wave_lens*Laser.k_prop) .* Mirror.mask;






