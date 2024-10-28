function SaveTempData(Params, Var, VarName)

N_DEG = Params.N_DEG;
gamma = Params.gamma;

file_name = ['.\Results\Temp\Temp_N_DEG_' strrep(strrep(num2str(N_DEG, '%1.1e'),'.','_'),'+','') '_gamma_'  strrep(strrep(num2str(gamma, '%1.1e'),'.','_'),'+','') '.mat'];

eval([ VarName ' = Var;' ]);
if (~exist(file_name, 'file'))
    save(file_name, 'Params');
end
eval(['save '  file_name ' -append ' VarName ';' ]);


