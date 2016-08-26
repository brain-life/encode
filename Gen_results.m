function [] = Gen_results(dataRootPath, dataOutputPath, subject, algorithm, parameter, fg_struct, connectome)

feFileName    = 'life_build_model'; 

dwiFile       = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
dwiFileRepeat = deblank(ls(fullfile(dataRootPath,subject,'diffusion_data','*b2000_aligned*.nii.gz')));
t1File        = deblank(fullfile(dataRootPath,subject,'anatomy',  'T1w_acpc_dc_restore_1p25.nii.gz'));

% Number of iterations for the optimization
Niter = 500;

L = 360;
tic

%% Build the model

fe = feConnectomeInit(dwiFile,fg_struct.fg,feFileName,[],dwiFileRepeat,t1File,L,[1,0],0); % We set dwiFileRepeat =  run 02
clear 'fg_struct'
disp(' ')
disp(['Time for model construction ','(L=',num2str(L),')=',num2str(toc),'secs']);

%% Fit the LiFE_SF model (Optimization)
tic
fe = feSet(fe,'fit',feFitModel(fe.life.M,feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));
disp(' ')

STC_run01.FittingTime = toc;
STC_run01.w = feGet(fe,'fiber weights');
STC_run01.subject = subject;
STC_run01.algorithm = algorithm;
STC_run01.parameter = parameter;
disp(['Time fitting LiFE=',num2str(STC_run01(1).FittingTime)]);
STC_run01.rmse = feGet(fe,'vox rmse');
STC_run01.rmsexv = feGetRep(fe,'vox rmse');
STC_run01.rrmse = feGetRep(fe,'vox rmse ratio');

disp('SAVING RESULTS...')
% save(fullfile(dataOutputPath,sprintf('Results_STC_run01_%s_%s_%s_conn%s.mat',subject,char(algorithm),char(parameter),connectome)), 'STC_run01','-v7.3')
% save(fullfile(dataOutputPath,sprintf('fe_structure_%s_%s_%s_%s_conn%s.mat',subject,'STC_run01',char(algorithm),char(parameter),connectome)), 'fe','-v7.3')


end







