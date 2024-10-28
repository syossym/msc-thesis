% Macro to find the resonance length of the cavity
% Do the phase average over some round trips

% Phase_shift = exp(sqrt(-1)*k_prop* Dist_reson);

Length.nb_propa_field = 200;        % Nb of light round trip to take into account

fprintf('\n------ Calculating the resonance length ------ \n')

Field.Circ = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Length.phase = zeros(1,Length.nb_propa_field,'double');
Grid.half_num = Grid.Num_point/2;


Field.Circ = Propa_mirror(Field.Start, Mirror.ITM_trans,i*ITM.t);

for q = 1:Length.nb_propa_field
    Length_tmp_phase = angle(Field.Circ(Grid.half_num,Grid.half_num));
    Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
    Field.Circ = Propa_mirror(Field.Circ,Mirror.ETM_cav,ETM.r);
    Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
    Field.Circ = Propa_mirror(Field.Circ, Mirror.ITM_cav,ITM.r);   
    Length.phase(q) = angle(Field.Circ(Grid.half_num,Grid.half_num))-Length_tmp_phase;
end
Length.phase = unwrap(Length.phase);
Phase_shift = exp(-i*mean(Length.phase));

fprintf('Finished ... \n')



