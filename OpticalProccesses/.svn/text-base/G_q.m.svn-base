%
% This function calculates the 2D form factor
%
function G = G_q(Xi12,k_vec,z_grid)

z_exp = -abs(repmat(z_grid,1,length(z_grid)) - repmat(z_grid',length(z_grid),1));

for (kk = 1:length(k_vec))
    z_exp_q = exp(z_exp.*k_vec(kk));
    integrand_G = Xi12.*z_exp_q;
    G(kk) = trapz(z_grid, trapz(z_grid,integrand_G,1), 2);
end
