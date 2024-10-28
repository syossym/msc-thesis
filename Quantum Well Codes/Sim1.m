clear all; close all; clc;

%% Exciton binding energy

x = 0.3;
% EX0_vec = [];
% lambda_vec = [];
% beta_vec = [];

%for (W = 5:5:150)
   W = 200;
   
   % calculate the wavefunctions
   [energies,wf,rts] = KP_WaveFunctions(x,W);
   %[beta,EX0,EX0_l] = ExcitonBindingEnergy(wf.wf_c(:,1), abs(wf.wf_v.lh{1}(:,1)));
   wf_hh_2 = [wf.z', abs(wf.wf_v.hh{1}(:,1))];
   wf_lh_1 = [wf.z', (wf.wf_v.lh{1}(:,1))];
   wf_e_1 = [wf.z', wf.wf_c(:,1)];
   
   figure(1);
   plot(wf_e_1(:,1), wf_e_1(:,2),'*-');
%    hold on; 
%    plot(wf_lh_1(:,1), wf_lh_1(:,2),'.-','Color','r');
   hold on; 
   plot(wf_hh_2(:,1), wf_hh_2(:,2),'o-','Color','g');
   
   pause;
   
   % delete all saved files
   
   dir_name = 'C:\Users\Yossi Michaeli\Documents\Thesis\Code\Quantum Well Codes\C Codes';
   delete(strcat(dir_name,'\*.r'));
   
   % save wave functions to files
   h_filename = strcat(dir_name,'\wf_h1.r');
   e_filename = strcat(dir_name,'\wf_e1.r');
   save(e_filename,'wf_e_1','-ascii');
   save(h_filename,'wf_hh_2','-ascii');
   %save(h_filename,'wf_hl_1','-ascii');
   
   % run the script
   
   open(strcat(dir_name,'\ebe.exe'))
   disp('Press any key to continue...');
   pause;
   
   % read the results
   M = dlmread(strcat(dir_name,'\EX0.r'));  

   % save results
   EX0_vec = [EX0_vec,M(1)];
   lambda_vec = [lambda_vec, M(2)];
   beta_vec = [beta_vec, M(3)];
   
%end