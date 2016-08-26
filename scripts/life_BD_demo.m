function [fh, fe] = life_BD_demo()
%% Example of initialization and fitting of the Big Data (BD) LiFE model (LiFE-BD)
% This demo function illustrates how to:
%  - A - Set up a LiFE-BD structure, identified as 'fe' (fascicle evaluation) in
%  the code below. This model contains a prediction of the diffusion
%  measurements in each white-matter voxel given by the fascicles contained
%  in a tractogrpahy solution, the connectome. Each fascicle makes a
%  prediction about the direction of diffusion in the set of voxels where
%  it travels through. The prediction is generated given the fascicle
%  orientation and position inside the voxel. Predictions from multiple
%  fascicles in each voxels are combined to generate a global connectome
%  prediciton for the diffusion signal in large sets of white matter
%  voxels.
%  - B - Fit the LiFE-BD model to compute the weights associated to each fascicle
%  in the connectome. Fascicles in the conenctome contribute differently to
%  predicting the diffusion signal in each voxel. First of all, fascicles
%  make predictions about the diffusion only in voxels where they travel.
%  Secondly, some fascicles have paths that produce better diffusion
%  predictions than others. We use a least-square method to find the
%  contribution of each fascicle to the diffusion signal in the white
%  matter voxels where the fascicles travels. A single weight is assigned
%  to each fascicle representing the global contribution of the fasicle to
%  the signal of all the voxels along its path - we call this
%  fascicle-global. Because multiple fascicles exist in several voxels the
%  set of fascicles weights and fascicles predicitons represents the
%  connectome-global prediction of the diffusion signal within the entire
%  set of white matter voxels. Estimating the fascicle weights allows for
%  evaluating the quality of the tractography solution. Eliminating
%  fascicles that do not contribute to predicting the diffusion signal
%  (they have assigned a zreo-weight). Finaly, the root-mean-squared error
%  (RMSE) of the model to the diffusion data - the model prediction error -
%  is used to evaluate the model prediction quality, compare different
%  tractography models and to perform statistical inference on the on
%  properties of the connectomes.
%  - C - Compare two different connectome models. This demo will show how to
%  compare two different conenctome models by using the diffusion
%  prediction error (the Root-Mean-Squared Error, RMSE). We report the
%  example of two conenctomes one generated using Constrained-spherical
%  deconvolution (CSD) and probabilistic tractography the other using a
%  tensor model and deterministic tractography
%  - Note - The example connectomes used for this demo comprise a portion
%  of the right occiptial lobe of an individual human brain. LiFE-BD utilizes
%  large-scale methods optimized for an efficient use of memory to solve the 
%  foward model. The software allows for solving connectomes spanning the 
%  entire white-matter of idnvidual brains. The size of the connectome on 
%  the test data set is small enought to allow for testing the code within 
%  a few minutes requiring only about 10GB of computer RAM and standard 
%  hardaware. This code has been tested with: 
%

%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Intialize a local matlab cluster if the parallel toolbox is available.
% This helps speeding up computations espacially for large conenctomes.

%feOpenLocalCluster;

%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'life_demo_anatomy_t1w_stanford.nii.gz');

%% (1) Evaluate the Probabilistic CSD-based connectome.
% We will analyze first the CSD-based probabilistic tractography
% connectome.
prob.tractography = 'Probabilistic';
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'LiFE-BD_build_model_demo_CSD_PROB';

%% (1.1) Initialize the LiFE-BD model structure, 'fe' in the code below. 
% This structure contains the forward model of diffusion based on the
% tractography solution. It also contains all the information necessary to
% compute model accuracry, and perform statistical tests. You can type
% help('feBuildModel') in the MatLab prompt for more information.

N = 360; % Discretization parameter

mycomputer = computer();
release = version('-release');
switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a','MACI64_2014b'}
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],0);
        otherwise
        sprintf('WARNING: currently LiFE is optimized for an efficient usage of memory \n using the Sparse Tucker Decomposition aproach (Caiafa&Pestilli, 2015) \n ONLY for Linux (MatlabR2015a) and MacOS (MatlabR2014b). \n If you have a different system or version you can still \n use the old version of LiFE (memory intensive). \n\n')
        sprintf('\n Starting building big matrix M in OLD LiFE...\n')
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],1);
end

%% (1.2) Fit the model. 
% Hereafter we fit the forward model of tracrography using a least-squared
% method. The information generated by fitting the model (fiber weights
% etc) is then installed in the LiFE-BD structure.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% (1.3) Extract the RMSE of the model on the fitted data set. 
% We now use the LiFE-BD structure and the fit to compute the error in each
% white-matter voxel spanned by the tractography model.
prob.rmse   = feGet(fe,'vox rmse');

%% (1.4) Extract the RMSE of the model on the second data set. 
% Here we show how to compute the cross-valdiated RMSE of the tractography
% model in each white-matter voxel. We store this information for later use
% and to save computer memory.
prob.rmsexv = feGetRep(fe,'vox rmse');

%% (1.5) Extract the Rrmse. 
% We show how to extract the ratio between the model prediction error
% (RMSE) and the test-retest reliability of the data.
prob.rrmse  = feGetRep(fe,'vox rmse ratio');

%% (1.6) Extract the fitted weights for the fascicles. 
% The following line shows how to extract the weight assigned to each
% fascicle in the connectome.
prob.w      = feGet(fe,'fiber weights');

%% (1.7) Plot a histogram of the RMSE. 
% We plot the histogram of  RMSE across white-mater voxels.
[fh(1), ~, ~] = plotHistRMSE(prob);

%% (1.8) Plot a histogram of the RMSE ratio.
% As a reminder the Rrmse is the ratio between data test-retest reliability
% and model error (the quality of the model fit).
[fh(2), ~] = plotHistRrmse(prob);

%% (1.9) Plot a histogram of the fitted fascicle weights. 
[fh(3), ~] = plotHistWeights(prob);
clear fe

switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a','MACI64_2014b'}
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],0);
        otherwise
        sprintf('WARNING: currently LiFE is optimized for an efficient usage of memory \n using the Sparse Tucker Decomposition aproach (Caiafa&Pestilli, 2015) \n ONLY for Linux (MatlabR2015a) and MacOS (MatlabR2014b). \n If you have a different system or version you can still \n use the old version of LiFE (memory intensive). \n\n')
        sprintf('\n Starting building big matrix M in OLD LiFE...\n')
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],1);
end


%% Extract the coordinates of the white-matter voxels
% We will use this later to compare probabilistic and deterministic models.
p.coords = feGet(fe,'roi coords');
clear fe

%% (2) Evaluate the Deterministic tensor-based connectome.
% We will now analyze the tensor-based Deterministic tractography
% connectome.
det.tractography = 'Deterministic';
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'life_demo_mrtrix_tensor_deterministic.mat');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'LiFE-BD_build_model_demo_TENSOR_DET';

%% (2.1) Initialize the LiFE-BD model structure, 'fe' in the code below. 
% This structure contains the forward model of diffusion based on the
% tractography solution. It also contains all the information necessary to
% compute model accuracry, and perform statistical tests. You can type
% help('feBuildModel') in the MatLab prompt for more information.
clear fe
switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a','MACI64_2014b'}
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],0);
        otherwise
        sprintf('WARNING: currently LiFE is optimized for an efficient usage of memory \n using the Sparse Tucker Decomposition aproach (Caiafa&Pestilli, 2015) \n ONLY for Linux (MatlabR2015a) and MacOS (MatlabR2014b). \n If you have a different system or version you can still \n use the old version of LiFE (memory intensive). \n\n')
        sprintf('\n Starting building big matrix M in OLD LiFE...\n')
        fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,N,[1,0],1);
end



%% (2.2) Fit the model. 
% Hereafter we fit the forward model of tracrography using a least-squared
% method. The information generated by fitting the model (fiber weights
% etc) is then installed in the LiFE-BD structure.
fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% (2.3) Extract the RMSE of the model on the fitted data set. 
% We now use the LiFE-BD structure and the fit to compute the error in each
% white-matter voxel spanned by the tractography model.
det.rmse   = feGet(fe,'vox rmse');

%% (2.4) Extract the RMSE of the model on the second data set. 
% Here we show how to compute the cross-valdiated RMSE of the tractography
% model in each white-matter voxel. We store this information for later use
% and to save computer memory.
det.rmsexv = feGetRep(fe,'vox rmse');

%% (2.5) Extract the Rrmse. 
% We show how to extract the ratio between the model prediction error
% (RMSE) and the test-retest reliability of the data.
det.rrmse  = feGetRep(fe,'vox rmse ratio');

%% (2.6) Extract the fitted weights for the fascicles. 
% The following line shows how to extract the weight assigned to each
% fascicle in the connectome.
det.w      = feGet(fe,'fiber weights');

%% (2.7) Plot a histogram of the RMSE. 
% We plot the histogram of  RMSE across white-mater voxels.
[fh(1), ~, ~] = plotHistRMSE(det);

%% (2.8) Plot a histogram of the RMSE ratio.
% As a reminder the Rrmse is the ratio between data test-retest reliability
% and model error (the quality of the model fit).
[fh(2), ~] = plotHistRrmse(det);

%% (2.9) Plot a histogram of the fitted fascicle weights. 
[fh(3), ~] = plotHistWeights(det);

%% Extract the coordinates of the white-matter voxels.
% We will use this later to compare probabilistic and deterministic models.
d.coords = feGet( fe, 'roi coords');
clear fe

%% (3) Compare the quality of fit of Probabilistic and Deterministic connectomes.
%% (3.1) Find the common coordinates between the two connectomes.
%
% The two tractography method might have passed through slightly different
% white-matter voxels. Here we find the voxels where both models passed. We
% will compare the error only in these common voxels. There are more
% coordinates in the Prob connectome, because the tracking fills up more
% White-matter.
%
% So, hereafter:
% - First we find the indices in the probabilistic connectome of the
% coordinate in the deterministic connectome. But there are some of the
% coordinates in the Deterministic conectome that are NOT in the
% Probabilistic connectome.
%
% - Second we find the indices in the Deterministic connectome of the
% subset of coordinates in the Probabilistic connectome found in the
% previous step.
%
% - Third we find the common voxels. These allow us to find the rmse for
% the same voxels.
fprintf('Finding common brain coordinates between P and D connectomes...\n')
prob.coordsIdx = ismember(p.coords,d.coords,'rows');
prob.coords    = p.coords(prob.coordsIdx,:);
det.coordsIdx  = ismember(d.coords,prob.coords,'rows');
det.coords     = d.coords(det.coordsIdx,:);
prob.rmse      = prob.rmse( prob.coordsIdx);
det.rmse       = det.rmse( det.coordsIdx);
clear p d

%% (3.2) Make a scatter plot of the RMSE of the two tractography models
fh(4) = scatterPlotRMSE(det,prob);

%% (3.3) Compute the strength-of-evidence (S) and the Earth Movers Distance.
% Compare the RMSE of the two models using the Stregth-of-evidence and the
% Earth Movers Distance.
se = feComputeEvidence(prob.rmse,det.rmse);

%% (3.4) Strength of evidence in favor of Probabilistic tractography. 
% Plot the distributions of resampled mean RMSE
% used to compute the strength of evidence (S).
fh(5) = distributionPlotStrengthOfEvidence(se);

%% (3.5) RMSE distributions for Probabilistic and Deterministic tractography. 
% Compare the distributions using the Earth Movers Distance.
% Plot the distributions of RMSE for the two models and report the Earth
% Movers Distance between the distributions.
fh(6) = distributionPlotEarthMoversDistance(se);

end

% ---------- Local Plot Functions ----------- %
function [fh, rmse, rmsexv] = plotHistRMSE(info)
% Make a plot of the RMSE:
rmse   = info.rmse;
rmsexv = info.rmsexv;

figName = sprintf('%s - RMSE',info.tractography);
fh = mrvNewGraphWin(figName);
[y,x] = hist(rmse,50);
plot(x,y,'k-');
hold on
[y,x] = hist(rmsexv,50);
plot(x,y,'r-');
set(gca,'tickdir','out','fontsize',16,'box','off');
title('Root-mean squared error distribution across voxels','fontsize',16);
ylabel('number of voxels','fontsize',16);
xlabel('rmse (scanner units)','fontsize',16);
legend({'RMSE fitted data set','RMSE cross-validated'},'fontsize',16);
end

function [fh, R] = plotHistRrmse(info)
% Make a plot of the RMSE Ratio:

R       = info.rrmse;
figName = sprintf('%s - RMSE RATIO',info.tractography);
fh      = mrvNewGraphWin(figName);
[y,x]   = hist(R,linspace(.5,4,50));
plot(x,y,'k-','linewidth',2);
hold on
plot([median(R) median(R)],[0 1200],'r-','linewidth',2);
plot([1 1],[0 1200],'k-');
set(gca,'tickdir','out','fontsize',16,'box','off');
title('Root-mean squared error ratio','fontsize',16);
ylabel('number of voxels','fontsize',16);
xlabel('R_{rmse}','fontsize',16);
legend({sprintf('Distribution of R_{rmse}'),sprintf('Median R_{rmse}')});
end

function [fh, w] = plotHistWeights(info)
% Make a plot of the weights:

w       = info.w;
figName = sprintf('%s - Distribution of fascicle weights',info.tractography);
fh      = mrvNewGraphWin(figName);
[y,x]   = hist(w( w > 0 ),logspace(-5,-.3,40));
semilogx(x,y,'k-','linewidth',2)
set(gca,'tickdir','out','fontsize',16,'box','off')
title( ...
    sprintf('Number of fascicles candidate connectome: %2.0f\nNumber of fascicles in optimized connetome: %2.0f' ...
    ,length(w),sum(w > 0)),'fontsize',16)
ylabel('Number of fascicles','fontsize',16)
xlabel('Fascicle weight','fontsize',16)
end

function fh = scatterPlotRMSE(det,prob)
figNameRmse = sprintf('prob_vs_det_rmse_common_voxels_map');
fh = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([det.rmse;prob.rmse]',{[10:1:70], [10:1:70]});
ymap = ymap./length(prob.rmse);
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
ylabel('Deterministic_{rmse}','fontsize',16)
xlabel('Probabilistic_{rmse}','fontsize',16)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ...
    'yTickLabel',round(1000*10.^[min(tck),...
    mean(tck), ...
    max(tck)])/1000, ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16,'visible','on')
end

function fh = distributionPlotStrengthOfEvidence(se)

y_e        = se.s.unlesioned_e;
ywo_e      = se.s.lesioned_e;
dprime     = se.s.mean;
std_dprime = se.s.std;
xhis       = se.s.unlesioned.xbins;
woxhis     = se.s.lesioned.xbins;

histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('Strength_of_Evidence_test_PROB_vs_DET_model_rmse_mean_HIST');
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1});
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 .2], ... 
        'xlim',[min(xhis) max(woxhis)], ...
        'xtick',[min(xhis) round(mean([xhis, woxhis])) max(woxhis)], ...
        'ytick',[0 .1 .2], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

title(sprintf('Strength of evidence:\n mean %2.3f - std %2.3f',dprime,std_dprime), ...
    'FontSize',16)
legend({'Probabilistic','Deterministic'})
end

function fh = distributionPlotEarthMoversDistance(se)

prob = se.nolesion;
det  = se.lesion;
em   = se.em;

histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('EMD_PROB_DET_model_rmse_mean_HIST');
fh = mrvNewGraphWin(figName);
plot(prob.xhist,prob.hist,'r-','color',histcolor{1},'linewidth',4);
hold on
plot(det.xhist,det.hist,'r-','color',histcolor{2},'linewidth',4); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 .12], ... 
        'xlim',[0 95], ...
        'xtick',[0 45 90], ...
        'ytick',[0 .06 .12], ...
        'fontsize',16)
ylabel('Proportion white-matter volume','fontsize',16)
xlabel('RMSE (raw MRI scanner units)','fontsize',16')
title(sprintf('Earth Movers Distance: %2.3f (raw scanner units)',em.mean),'FontSize',16)
legend({'Probabilistic','Deterministic'})
end

