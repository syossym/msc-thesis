function Field_built = Build_Field_Cavity(Delta_L)

global Grid;
global Length;
global Field;
global Laser;

Field_built = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
%Field_built3 = zeros(Grid.Num_point,Grid.Num_point,Length.nb_propa_field,'double');
for q = 1:Length.nb_propa_field
    Field_built = Field_built + Field.propa(:,:,q) * exp(i*Laser.k_prop* Delta_L*q);
end

% tmp_shift(1,1,:) = exp(i*Laser.k_prop* 2 *Delta_L*(1:Length.nb_propa_field));
% Field_shift = repmat(tmp_shift,[Grid.Num_point Grid.Num_point 1]); 
% Field_built = dot(Field.propa,Field_shift,3);

% for q = 1:Length.nb_propa_field
%      Field_built3(:,:,q) = Field.propa(:,:,q) .* exp(i*Laser.k_prop* Delta_L*q);
%      Field_built = sum(Field_built3,3);
%  end

return