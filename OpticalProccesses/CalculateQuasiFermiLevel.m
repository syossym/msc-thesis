function Eqf = CalculateQuasiFermiLevel(type, E, N, T, m)

global Consts;

% Build input files - levels and concentration
E_file = [1 E/Consts.e_0/1e-3];
E_fid = fopen(['.\Bin\E' type '.r'], 'wt');
fprintf(E_fid, '%2i %E', E_file);
fclose(E_fid);

N_file = [1,N];
N_fid = fopen('.\Bin\N.r', 'wt');
fprintf(N_fid, '%2i %E', N_file);
fclose(N_fid);

% Update the potential profile
delete('.\Bin\v.r');
copyfile(['.\Bin\Output\v_' type '.r'], '.\Bin\v.r');

% Run the calcuation program
delete('.\Bin\Ef.r');
if (nargin == 5)
    RunExternalCode('.\Bin\sbp.exe', [' -p ' type ' -f true -T ' num2str(T) ' -m ' num2str(m)]);
else
    RunExternalCode('.\Bin\sbp.exe', [' -p ' type ' -f true -T ' num2str(T)]);
end
copyfile('.\Bin\Ef.r', ['.\Bin\Output\Ef_' type '.r']);
copyfile('.\Bin\FD*.r', '.\Bin\Output\');

% Get the result
Eqf_file = load( ['.\Bin\Output\Ef_' type '.r']);
Eqf = Eqf_file(1,2)*(Consts.e_0*1e-3);