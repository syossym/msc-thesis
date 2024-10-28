% Propagate the field inside the cavity for an arbituary microscopique cavity length 

Field.Circ = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.propa = complex(zeros(Grid.Num_point,Grid.Num_point,Length.nb_propa_field,'double'));

fprintf(' Calculating the field in the cavity ')

Field.propa(:,:,1) = Field.Start;
Field.Circ = Field.Start;
fprintf('\n Be patient ... calculating ...       ')
 for q = 2:(Length.nb_propa_field)
         
 	Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
 	Field.Circ = Propa_mirror(Field.Circ, Mirror.ETM_cav,ETM.r);
 	Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
 	Field.Circ = Propa_mirror(Field.Circ, Mirror.ITM_cav,ITM.r);
    Field.propa(:,:,q) = Field.Circ;
    
    if (rem(q,(Length.nb_propa_field)/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*q/(Length.nb_propa_field))       
    end     
 end

 fprintf('\b\b\b\b\b\b\b\b  Finished! ')  
 fprintf('\n')
