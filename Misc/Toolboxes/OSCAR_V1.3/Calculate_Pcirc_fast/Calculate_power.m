function power_temp = Calculate_power(wave_name)

global Grid;

power_temp = sum(sum(abs(wave_name).^2)) * Grid.step^2;



