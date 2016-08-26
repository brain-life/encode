function [fh, fe] = demo3_virtual_lesion()
%% Example of Virtual Lesion computation using the multidimensional LiFE model
% This demo function illustrates how to:
%  
%

%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa
%  (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Check if vistasoft is installed
check = which('vistaRootPath');
if isempty(check)
    disp('ERROR: vistasoft package not installed.')
    disp('Please, download it from https://github.com/vistalab/vistasoft and add it to the Matlab Path')
    return
end

% Check if mba is installed
check = which('mbaComputeFibersOutliers');
if isempty(check)
    disp('ERROR: mba package not installed.')
    disp('Please, download it from https://github.com/francopestilli/mba')
    return
end

%% Compute Virtual Lesion (VL) from fe structures.

% Change to datasets directory
[folder, name, ext] = fileparts(which(mfilename));
cd(folder);
cd ../demo_datasets/fe_structures/

% Choose major tract for computing VL
clc
prompt = 'Please select a major tract number (1 to 20): \n1-2: Anterior thalamic radiation (ATR) \n3-4: Cortico Spinal Tract (CST) \n5-6: Cingulum (cingulate gyrus) (Cing) \n7-8: Cingulum (hippocampus)  (Hipp) \n9-10: Forceps minor/major \n11-12: Inferior fronto-occipital fasciculus (InFOF) \n13-14: Inferior longitudinal fasciculus (InLF) \n15-16: Superior longitudinal fasciculus (SuLF) \n17-18: Uncinate fasciculus (UF) \n19:20: Superior longitudinal fasciculus (temporal part) (Temp) \n\n';
tract = input(prompt);

%% Read STN subject PROB results
disp('loading fe_structures for FP subject in STN dataset ...')
% load fe_structure
load('fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_PROB_lmax10_connNUM01.mat')
ind = find(~isnan(fe.life.fit.weights)&fe.life.fit.weights>0);

% load tract classification
load('fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_500000_SD_PROB_lmax10_connNUM01_TRACTS.mat')

tract_name = char(classification.names(tract));
disp(char(strcat('Computing virtual lesion of ', tract_name,' for subject in STN dataset (PROB)')));

ind_tracts1 = find(classification.index==tract); % indices to fascicles in the selected major tract

% We clean all major tracs
sprintf('\n Cleaning %s ...',tract_name)
[fg_tract, keep_tract] = mbaComputeFibersOutliers(fascicles(tract),3,3);
    
ind_tracts1 = ind_tracts1(keep_tract);

ind_nnz = find(fe.life.fit.weights);
ind1 = ind_nnz(ind_tracts1);
if isempty(ind1)
    rmse_wVL = [];
    rmse_woVL = [];
    nFib_tract = 0;
    nFib_PN = 0;
    nVoxels = 0;
else
    [rmse_wVL, rmse_woVL,nFib_tract, nFib_PN, nVoxels]= feComputeVirtualLesion_norm(fe,ind1);
end

if nFib_tract ==0
    se = [];
else
    se = feComputeEvidence_norm(rmse_woVL,rmse_wVL);
end


%% Strength of evidence in favor of Probabilistic tractography. 
% Plot the distributions of resampled mean RMSE
% used to compute the strength of evidence (S).
fh(1) = distributionPlotStrengthOfEvidence(se);

%% RMSE distributions for Probabilistic and Deterministic tractography. 
% Compare the distributions using the Earth Movers Distance.
% Plot the distributions of RMSE for the two models and report the Earth
% Movers Distance between the distributions.
fh(2) = distributionPlotEarthMoversDistance(se);


%% Generate visualization of the major tract and its path neighborhood
colors     = {[.1 .25 .65],[.75 .25 .1]};
viewCoords = [90,0];
slice      = [-1 0 0];
proportion_to_show = .05;
threshold_length = 10;


ind_nnz = find(fe.life.fit.weights);
ind1    = ind_nnz(ind_tracts1);
ind_tracts2 = feGet(fe,'Path Neighborhood',ind1);
fg          = feGet(fe,'fibers img');
fg_pathn    = fgExtract(fg,ind_tracts2,'keep');
clear fg
roicoords = feGet(fe,'coords from fibers',ind1);

fg_pathn.fibers = mbaFiberSplitLoops(fg_pathn.fibers);
xform     = feGet(fe,'img2acpcxform');
fg_pathn  = dtiXformFiberCoords(fg_pathn,xform,'acpc');

% Prepare the plot of tract and neighborhood
fg{1}    = fg_tract;
fg_pathn = rmfield(fg_pathn,'coordspace');
fg_pnplot = fg_pathn;

c = 1;
for ii = 1:length(fg_pathn.fibers)
    if length(fg_pathn.fibers{ii}) > threshold_length
        fibers{c} = fg_pathn.fibers{ii};
        c = c + 1;
    end
end
fg_pathn.fibers = fibers; clear fibers
fibs_indx = randsample(1:length(fg_pathn.fibers),round(length(fg_pathn.fibers)*proportion_to_show));
fg_pnplot.fibers = fg_pnplot.fibers(fibs_indx);
fg{2}    = fg_pnplot;

% plot tract and PN
fig_name      = char(strcat(tract_name,'+ PN (only ',num2str(proportion_to_show*100),'% of PN fascicles)'));
[fh(3), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1 2]);
%feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
%close all; drawnow

% plot tract
fig_name      = char(strcat(tract_name));
[fh(4), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1]);
%feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
%close all; drawnow

% plot PN
fig_name      = char(strcat('PN (only ',num2str(proportion_to_show*100),'% of PN fascicles)') );
[fh(5), ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [2]);
%feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
%close all; drawnow



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
