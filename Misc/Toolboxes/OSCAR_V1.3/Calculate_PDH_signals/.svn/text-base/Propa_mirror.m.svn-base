function Output = Propa_mirror(Wave_field, Wave_mirror, reflec)

global Mirror
global Laser;

Output = Wave_field .* exp(-i * Wave_mirror*Laser.k_prop) * reflec .* Mirror.mask;

