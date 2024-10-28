function Field_built = Build_Field_Cavity(Delta_L)

global Grid;
global Length;
global Field;
global Laser;

Field_built = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));

for q = 1:Length.nb_propa_field
    Field_built = Field_built + Field.propa(:,:,q) * exp(i*Laser.k_prop* Delta_L*q);
end

return