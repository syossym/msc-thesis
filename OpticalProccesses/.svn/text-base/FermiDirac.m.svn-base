
function fd = FermiDirac(E_f, E, T)

global Consts;

fd = 1./(exp((E - E_f*ones(1,length(E)))/(Consts.k_B*T)) + 1);
