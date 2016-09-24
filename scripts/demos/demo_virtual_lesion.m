function [fh, fe] = demo_virtual_lesion()
% Example of Virtual Lesion computation using the multidimensional encoding
% model and the LiFE method.
% 
% This demo function illustrates how to perfomr a virtual lesion by using
% the multidimensional connectome encoding framework.
%
% The demo reproduces some of the results initially published in Pestilli
% et al., Nature Methods 2014 and replicated in Caiafa and Pestilli
% forthcoming.
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa
%  (CONICET)
%
%  email: frakkopesto@gmail.com and ccaiafa@gmail.com

%% (0) Check matlab dependencies and path settings.
if ~exist('vistaRootPath.m','file');
    disp('Vistasoft package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/vistalab/vistasoft');
end
if ~exist('mbaComputeFibersOutliers','file')
    disp('ERROR: mba package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/francopestilli/mba')
end
if ~exist('feDemoDataPath.m','file');
    disp('ERROR: demo dataset either not installed or not on matlab path.')
    error('Please, download it from http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz')
end

%% (1) Compute Virtual Lesion (VL) from fe structures.
%
% The demo data set provides a series of precomputed fascicles. These
% fascicles were segmented given known anatomical atlases of thehuman white
% matter (Mori, Susumu, et al. MRI atlas of human white matter. Elsevier,
% 2005.) using the AFQ toolbox (https://github.com/jyeatman/AFQ).
%
% The demo will ask to select one out of 20 major tracts. We will then use
% that tract to perform a virtual leasion using the LiFE toolbox. The
% virtual lesion operation will inform us of the statistical evidence for
% the tract given the tractogprahy solution and the data set provided.
%

% We load one precomputed LiFE structure (FE structure)
%
% The structure we load is provided as part of the Demo Data set.
%
% It was generated using one subject from the STN diffusion-weighted dataset 
% the MRTRIX toolbox, probabilistic tracking based of the CSD model (Lmax=10).
disp('loading fe_structures for FP subject in STN dataset ...')
feFileName = fullfile(feDemoDataPath('STN','sub-FP','fe_structures'), ...
'fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_PROB_lmax10_connNUM01.mat');
load(feFileName)

% All fascicles in a full brain connectome have been encoded into a three
% dimensional array, tensor Phi. Phi is available inside the FE sctructure
% we just loaded. Each dimension of Phi encodes different properties of the
% connectome. Mode 1 the orientation of the connectome fascicles. Mode 2
% the spatial location of each node/fascile. Mode 3 the identify of each
% fascicles, also called fascicle indices.
% 
% A white matter tract is defined as a set of fascicles. To identify a
% tract within Phi we need to select a group of fascicles within Phi along Mode 3.
%
% First we find the indices of all fascicles in the connectome that have a
% non-zero weight associated ('ind nzw'). These are fascicles that contributed a
% reliable amount in predicting the diffusion signal. A simple call of the
% hub function feGet.m helps with this.
ind_nzw = feGet(fe,'ind nzw');

% After that we load a precomputed tract segmentation (classification). A
% segmentaion of tracts was performed on the connectome. Edges of the
% connectome were classified as being part of one of twenty major human
% white matter tracts. For example the arcuate fasciculus, or the
% cortico-spinal tract.
% 
FileName = fullfile(feDemoDataPath('STN','sub-FP','tracts_classification'), ...
'fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_500000_SD_PROB_lmax10_connNUM01_TRACTS.mat');
load(FileName) % Load tract classification file from disk.

% Pick a tract and return indices of tract-fascicles in the encoded 
% connectome (Phi tensor).
[ind_fascicles_in_tract, ~, ~] = demo_local_choose_tract(fe,classification, fascicles, ind_nzw);

% Remove the fascicles of the tract to be lesioned (actually perform the
% lesion) from the encoded connectome. 
%
% Compute the root-mean-squared error of the model with and without tract
% in the white-matter voxels comprised by the tract. 
[rmse_wVL, rmse_woVL]= feComputeVirtualLesion_norm(fe,ind_fascicles_in_tract);

% Compute the statistical strenght of evidence for the tract. 
%
% This means measure the change in root-mean-squared error in predicting
% the measured demeaned diffusion-weighted MRI signal due to the tract
% lesion.
se = feComputeEvidence_norm(rmse_woVL,rmse_wVL);

% Plot the Strength of evidence. 
%
% Plot the distributions of resampled mean RMSE
% used to compute the strength of evidence (S).
fh(1) = distributionPlotStrengthOfEvidence(se);

% Plot the two RMSE distributions with and without lesion.
%
% Compare the distributions using the Earth Movers Distance. Plot the
% distributions of RMSE for the two models and report the Earth Movers
% Distance between the distributions.
fh(2) = distributionPlotEarthMoversDistance(se);

%% (3) Plot the anatomy of the tract and its path-neighborhood.
%
% Below we show how to extract a path neighborhood of a tract and the tract
% out of the FE structure, given only the indices of the fascicles of the
% tract in the condidate whole-brain tractography.

% Get full candidate connectome.
fg_whole_brain = feGet(fe,'fibers acpc');

% Get tract fascicles.
fg_tract       = fgExtract(fg_whole_brain,ind_fascicles_in_tract,'keep'); 

% Get path-neightborhood of the tract using the connectome encoding, via
% feGet.m
ind_pathneighborhood = feGet(fe,'Path Neighborhood',ind_fascicles_in_tract);
fg_pathn             = fgExtract(fg_whole_brain,ind_pathneighborhood,'keep');
clear fg_whole_brain

% Plot the anatomy of tract and neighborhood.
fh = [fh, demo_local_plot_anatomy(fg_tract, fg_pathn)];

end

%
% -- local helper functions -- %
%
function fh = demo_local_plot_anatomy(fg_tract,fg_pathn)
% - 
% Visualize the major tract and its path neighborhood
%
colors     = {[.1 .25 .65],[.75 .25 .1]};
viewCoords = [90,0];
proportion_to_show = .05;
threshold_length   = 10;

% Prepare the plot of tract and neighborhood
fg{1} = fg_tract; 

% We split the fibers of the path-neighborhood that enter and exit the
% tract ROI multiple times into separate fascicles. This is convenient for
% visualizing the fascicles.
fg_pathn.fibers = mbaFiberSplitLoops(fg_pathn.fibers);
fg_pathn = rmfield(fg_pathn,'coordspace');
c = 1;
for ii = 1:length(fg_pathn.fibers)
    if length(fg_pathn.fibers{ii}) > threshold_length
        fibers{c} = fg_pathn.fibers{ii};
        c = c + 1;
    end
end
fg_pathn.fibers = fibers; clear fibers
% Pick a percentage of fascicles to display (the PN can be too dense for visualization).
fibs_indx = randsample(1:length(fg_pathn.fibers), ...
            round(length(fg_pathn.fibers)*proportion_to_show));
fg{2}.fibers = fg_pathn.fibers(fibs_indx);

% plot tractPath-Neighborhood
fig_name      = char(strcat('Tract + PN (',num2str(proportion_to_show*100),'% of PN fascicles)'));
[fh(1), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1 2]);

dootherplots = true;
if dootherplots
% plot tract
fig_name      = 'Tract only';
[fh(2), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1]);

% plot PN
fig_name      = char(strcat('PN (only ',num2str(proportion_to_show*100),'% of PN fascicles)') );
[fh(3), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [2]);
end

% Change plot view.
view(0,90)

end



function [ind_good_fascicles_in_tract, tract_num, tract_name] = demo_local_choose_tract(fe,classification, fascicles, ind_nzw)
% 
% Local function to select a tract from a series of segemnted tracts.
%
% Aask to choose a major tract as example.
% Perform some preprocessing to return the appropriate tract indices that
% can be used for the virtual lesion process.
%
prompt = 'Please select a major tract number (1 to 20): \n1-2: Anterior thalamic radiation (ATR) \n3-4: Cortico Spinal Tract (CST) \n5-6: Cingulum (cingulate gyrus) (Cing) \n7-8: Cingulum (hippocampus)  (Hipp) \n9-10: Forceps minor/major \n11-12: Inferior fronto-occipital fasciculus (InFOF) \n13-14: Inferior longitudinal fasciculus (InLF) \n15-16: Superior longitudinal fasciculus (SuLF) \n17-18: Uncinate fasciculus (UF) \n19:20: Superior longitudinal fasciculus (temporal part) (Temp) \n\n';
tract_num = input(prompt);
tract_name  = char(classification.names(tract_num));
fprintf('[%s] Extracting tract %s from Encoding model... \n', mfilename, tract_name);

% First, given a selected tract index (tract_num), we find the indices of
% the fascicles in the cadidate connectome (all fascicles returned by
% tractography).
ind_tract_fascicles = find(classification.index == tract_num); 

% We clean all major tracs. This means that we remove anatomical outliers, namely fascicles too far from the mean anatomical location and lenght of the  This is because some of the initial
% segmentation performed by AFQ can return tracts that are too far away
% from the expected tract path. To overcome this limitation we accept
% tracts that are (1) close by the mean tract path coordinates, (2) not
% tool long or too short from the average length of the tract fascicles.
[~, fascicles_to_keep] = mbaComputeFibersOutliers(fascicles(tract_num),3,3);
fascicles_to_keep      = ind_tract_fascicles(fascicles_to_keep);

% Now we use the indices of the fascicles in the tract (ind_tracts1) and the indices of the
% non-zero weight fibers (ind_nzw) to identify the subset of tract 1 that
% is supported by the data.
ind_good_fascicles_in_tract   = ind_nzw(fascicles_to_keep);

end

function fh = distributionPlotStrengthOfEvidence(se)
%
% Make a plot of the Strengh-of-evidence result.
%
y_e        = se.s.unlesioned_e;
ywo_e      = se.s.lesioned_e;
dprime     = se.s.mean;
std_dprime = se.s.std;
xhis       = se.s.unlesioned.xbins;
woxhis     = se.s.lesioned.xbins;

histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('Strength_of_Evidence_test_without_VL_vs_with_VL');
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1});
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 1], ... 
        'xlim',[0 0.2], ...
        'xtick',[0 0.1 0.2], ...
        'ytick',[0 .5 1], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

title(sprintf('Strength of evidence:\n mean %2.3f - std %2.3f',dprime,std_dprime), ...
    'FontSize',16)
legend({'without VL','with VL'})
end

function fh = distributionPlotEarthMoversDistance(se)
%
% Make a plot fo the Earth Mover's distance result.
%
nolesion = se.nolesion;
lesion  = se.lesion;
em   = se.em;

histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];
figName = sprintf('EMD_without_VL_vs_with_VL');
fh = mrvNewGraphWin(figName);
plot(nolesion.xhist,nolesion.hist,'r-','color',histcolor{1},'linewidth',4);
hold on
plot(lesion.xhist,lesion.hist,'r-','color',histcolor{2},'linewidth',4); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ticklen',[.025 .05], ...
        'ylim',[0 .12], ... 
        'xlim',[0 0.2], ...
        'xtick',[0 0.1 0.2], ...
        'ytick',[0 .5 1], ...
        'fontsize',16)
ylabel('Proportion white-matter volume','fontsize',16)
xlabel('RMSE','fontsize',16')
title(sprintf('Earth Movers Distance: %2.3f',em.mean),'FontSize',16)
legend({'without VL','with VL'})
end

% Local functions to plot the tracts
function [fig_h, light_h] = plotFasciclesNoAnat(fascicles, color, viewCoords, fig_name,tracts_to_clean)
fig_h = figure('name',fig_name,'color','k');
hold on
set(gca,'visible','off','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78],'Color','w')
for iFas  = 1:length(tracts_to_clean)
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
%set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
set(gcf,'Color',[1 1 1])
drawnow

end
