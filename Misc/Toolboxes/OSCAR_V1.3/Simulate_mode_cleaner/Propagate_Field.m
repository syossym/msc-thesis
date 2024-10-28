% Propagate the field inside the cavity for an arbituary microscopique cavity length 

Field.Circ = complex(zeros(Grid.Num_point,Grid.Num_point,'double'));
Field.propa = complex(zeros(Grid.Num_point,Grid.Num_point,Length.nb_propa_field,'double'));

fprintf(' Calculating the field in the cavity ')

Field.Circ = Propa_mirror(Field.Start,Mirror.M1_trans,i*M1.t,Grid_beta.X);

Field.propa(:,:,1) = Field.Circ;

fprintf('\n Be patient ... calculating ...       ')
for q = 2:(Length.nb_propa_field)
  
     Field.Circ = Make_propagation(Field.Circ,MC.Length_small);
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M2_cav,M2.r,Grid_beta.X);
     Field.Circ = Make_propagation(Field.Circ,MC.Length_side);
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M3_cav,M3.r,Grid_alpha.X);
     Field.Circ = Make_propagation(Field.Circ,MC.Length_side);
     Field.Circ = Propa_mirror(Field.Circ,Mirror.M1_cav,M1.r,Grid_beta.X);

     Field.propa(:,:,q) = Field.Circ;
     
%  	Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
%  	Field.Circ = Propa_mirror(Field.Circ, Mirror.ETM_cav,ETM.r);
%  	Field.Circ = Make_propagation(Field.Circ,Mat_propagation);
%  	Field.Circ = Propa_mirror(Field.Circ, Mirror.ITM_cav,ITM.r);
%     Field.propa(:,:,q) = Field.Circ;
    
    if (rem(q,(Length.nb_propa_field)/100) == 0)
        fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*q/(Length.nb_propa_field))       
    end     
 end

 fprintf('\b\b\b\b\b\b\b\b  Finished! ')  
 fprintf('\n')
