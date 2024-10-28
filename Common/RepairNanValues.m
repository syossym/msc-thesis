function repaired = RepairNanValues(in_vec)

nan_indices = find(isnan(in_vec));

for (ii=1:length(nan_indices))
    if (nan_indices(ii) == 1)
        in_vec(nan_indices(ii)) = in_vec(nan_indices(ii)+1);
    elseif (nan_indices(ii) == length(in_vec))
        in_vec(nan_indices(ii)) = in_vec(nan_indices(ii)-1);  
    else
        if (in_vec(nan_indices(ii)+1)>in_vec(nan_indices(ii)-1))
            in_vec(nan_indices(ii)) = in_vec(nan_indices(ii)-1) + 0.5*(in_vec(nan_indices(ii)+1)-in_vec(nan_indices(ii)-1)); 
        else
            in_vec(nan_indices(ii)) = in_vec(nan_indices(ii)-1) - 0.5*(in_vec(nan_indices(ii)-1)-in_vec(nan_indices(ii)+1)); 
        end
    end
end

repaired = in_vec;