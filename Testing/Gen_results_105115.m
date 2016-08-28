function [] = Gen_results_105115()
%% Script to test on one subject HCP3T
clear
clc
%diary on

vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
addpath(genpath(vista_soft_path));

% Define path to the NEW LiFE
new_LiFE_path = '/N/dc2/projects/lifebid/code/ccaiafa/life/';
addpath(genpath(new_LiFE_path));

dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/';

subject = '105115';
conn = 'NUM01'; % 

param = 'lmax10'; % {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8', ''}
alg = 'SD_PROB'; % {'SD_PROB', 'SD_STREAM','tensor'}

fgFileName    = deblank(ls(fullfile(dataRootPath,subject,'fibers_new', strcat('*b2000*',char(param),'*',char(alg),'*',conn,'*','500000.tck'))));

fg = dtiImportFibersMrtrix(fgFileName);
fg_struct.fg = fg;

disp(['Gen results ', subject, '  connectome' ,conn,', ',alg,' ',param,'...'])
[fe_test, STC_run01_test] = Gen_results(dataRootPath, subject, alg, param, fg_struct, conn);

test.rmse   = feGet(fe_test,'vox rmse');
test.w      = feGet(fe_test,'fiber weights');
test.coords = feGet(fe_test,'roi coords');

disp(['Supported fibers (Orig)=',num2str(nnz(test.w))])
disp(['Global rmse (Orig)=',num2str(mean(test.rmse(:)))])


% load results computed for the paper
load(fullfile(dataOutputPath,sprintf('Results_STC_run01_%s_%s_%s_conn%s.mat',subject,char(alg),char(param),conn)))
load(fullfile(dataOutputPath,sprintf('fe_structure_%s_%s_%s_%s_conn%s.mat',subject,'STC_run01',char(alg),char(param),conn)))

orig.rmse   = feGet(fe,'vox rmse');
orig.w      = feGet(fe,'fiber weights');
orig.coords = feGet(fe,'roi coords');

disp(['Supported fibers (Orig)=',num2str(nnz(orig.w))])
disp(['Global rmse (Orig)=',num2str(mean(orig.rmse(:)))])

fprintf('Finding common brain coordinates between Test and Orig connectomes...\n')
test.coordsIdx = ismember(test.coords,orig.coords,'rows');
test.coords    = test.coords(test.coordsIdx,:);
orig.coordsIdx  = ismember(orig.coords,test.coords,'rows');
orig.coords     = orig.coords(orig.coordsIdx,:);
test.rmse      = test.rmse(test.coordsIdx);
orig.rmse       = orig.rmse(orig.coordsIdx);

fh = scatterPlotRMSE(test,orig);

rmpath(genpath(new_LiFE_path));
rmpath(genpath(vista_soft_path));

end


function [fe, STC_run01] = Gen_results(dataRootPath, subject, algorithm, parameter, fg_struct, connectome)

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

end


function fh = scatterPlotRMSE(orig,test)
figNameRmse = sprintf('orig_vs_test_rmse_common_voxels_map');
fh = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([orig.rmse;test.rmse]',{[10:10:1000], [10:10:1000]});
ymap = ymap./length(test.rmse);
sh   = imagesc(flipud(log10(ymap)));
cm   = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',16,'visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('orig_{rmse}','fontsize',16)
xlabel('test_{rmse}','fontsize',16)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16,'visible','on')
end










