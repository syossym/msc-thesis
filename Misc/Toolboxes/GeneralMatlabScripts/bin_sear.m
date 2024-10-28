function IndVal=bin_sear(X,Val);
%--------------------------------------------------------------------
% bin_sear function         binary search on a vector.
% input  : - sorted vector (ascending).
%	   - value to search.
% output : - index of closest value.
% Tested : Matlab 4.2
%          Matlab 7.0
%     By : Eran O. Ofek           September 1994
%                                      June 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
N      = length(X);

if (N==1),
   IndVal = 1;
else
   Ind1   = 1;
   Ind2   = N;
   IndM   = floor(0.5.*N);
   Y1     = X(Ind1);
   Y2     = X(Ind2);
   Ym     = X(IndM);

   Found  = 0;
   while (Found==0),
      %[[Ind1, IndM, Ind2]-499000, [Y1,   Ym,   Y2]]
      if (Val>Ym),
         Ind1 = IndM;
         %Ind2 = Ind2;
         Y1   = X(Ind1);
         %Y2   = X(Ind2);

         if ((Ind2-Ind1)>=2),
            IndM = floor(floor(0.5.*(Ind2+Ind1)));
         else
            Found = 1;
            if (abs(Val-Y1)<abs(Val-Y2)),
               IndVal = Ind1;
            else
               IndVal = Ind2;
            end
         end

         Ym   = X(IndM);

      elseif (Val<Ym),
         Ind2 = IndM;
         %Ind1 = Ind1;
         %Y1   = X(Ind1);
         Y2   = X(Ind2);

         if ((Ind2-Ind1)>=2),
            IndM = floor(floor(0.5.*(Ind2+Ind1)));
         else
            Found = 1;
            if (abs(Val-Y1)<abs(Val-Y2)),
               IndVal = Ind1;
            else
               IndVal = Ind2;
            end
         end

         Ym   = X(IndM);

      else
         Found  = 1;
         IndVal = IndM;
      end
   end
end


