function dirName = CreateSimulationResultsDir(parentDirPath, L_z, T, x, cb_len, vb_len, gamma)

dirName = ['Single_QW_FCT_HF_' num2str(L_z/1e-10) 'A_' num2str(T) 'K_x_' strrep(num2str(x), '.', '_') '_' num2str(cb_len) 'e_' num2str(vb_len) 'h_gamma_' num2str(gamma,'%1.0e')];
[status,message,messageid] = mkdir(parentDirPath, dirName);