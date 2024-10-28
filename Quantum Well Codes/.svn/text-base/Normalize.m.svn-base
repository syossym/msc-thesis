function f_norm = Normalize(f)

f_norm = f;
[f_max, mm] = max(f);
if (mm > length(f)*0.75) 
   f_norm(ceil(length(f)*0.75) : end)= 0;
end