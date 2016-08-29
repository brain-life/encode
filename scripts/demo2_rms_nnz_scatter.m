function [fh] = demo2_rms_nnz_scatter()
%% Example of comparing LiFE optimized connectomes
% This demo function reads LiFE multidimensional model from disk and
% display each subject in a scatter plot as function of the achieved root
% mean squated errors (rmse) and the obtained number of fascicles
% (Fascicles number). We show the results for 3 subjects, one belongin to
% each of the dataset used in the paper, i.e. STN (Stanford), HCP3T (Human 
% Connectome Project 3 Tesla) and HCP7T (Human Connectome Project 7 Tesla).
% Two examples of tratography algorithm are shown here: Probabilistic with
% Lmax=10 and Deterministic with Lmax=10.
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

fh = figure('name','combined scatter mean ±sem across repeats','color','w');
set(fh,'Position',[0,0,800,600]);
hold on

%% Read rmse (root mean squared error) and number of fibers nnz (number of non zero coeffs) from fe structures.

% Change to datasets directory
[folder, name, ext] = fileparts(which(mfilename));
cd(folder);
cd /N/dc2/projects/lifebid/code/ccaiafa/demo_datasets/fe_structures

%% Read STN subject PROB results
disp('loading fe_structures for FP subject in STN dataset ...')
load('fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_PROB_lmax10_connNUM01.mat')
sbj = retrieve_results(fe,'PROB','STN');
% plot point
Gen_plot(sbj,'medium')

%% Read STN subject DET results
load('fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_STREAM_lmax10_connNUM01.mat')
sbj = retrieve_results(fe,'DET', 'STN');
% plot point
Gen_plot(sbj,'medium')

%% Read HCP3T subject PROB results
disp('loading fe_structures for105115 subject in HCP3T dataset ...')
load('fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat')
sbj = retrieve_results(fe,'PROB', 'HCP3T');
% plot point
Gen_plot(sbj,'cold')

%% Read HCP3T subject DET results
load('fe_structure_105115_STC_run01_SD_STREAM_lmax10_connNUM01.mat')
sbj = retrieve_results(fe,'DET', 'HCP3T');
% plot point
Gen_plot(sbj,'cold')

%% Read HCP7T subject PROB results
disp('loading fe_structures for 108323 subject in HCP7T dataset ...')
load('fe_structure_108323_STC_run01_SD_PROB_lmax8_connNUM01.mat')
sbj = retrieve_results(fe,'PROB', 'HCP7T');
% plot point
Gen_plot(sbj,'hot')

%% Read HCP7T subject DET results
load('fe_structure_108323_STC_run01_SD_STREAM_lmax8_connNUM01.mat')
sbj = retrieve_results(fe,'DET', 'HCP7T');
% plot point
Gen_plot(sbj,'hot')

% Format figure
set(gca,'tickdir','out', 'ticklen',[0.025 0.025], ...
         'box','off','ytick',[2 9 16].*10^4, 'xtick', [0.04 0.07 0.1], ...
         'ylim',[2 16].*10^4, 'xlim', [0.04 0.1],'fontsize',20)
axis square
ylabel('Fascicles number','fontsize',20)
xlabel('Connectome error (r.m.s.)','fontsize',20)
legend1 = legend('show');
set(legend1,...
    'Position',[0.747708353369186 0.681944453087118 0.234999994523823 0.30083332469066]);


end

function [] = Gen_plot(sbj,color_type)

c = getNiceColors(color_type);

%% scatter plot
a = 0.5;

switch sbj.alg
    case 'PROB'
        plot(sbj.rmse, sbj.nnz,'o','markerfacecolor',c(1,:),'markeredgecolor','k','linewidth',a,'markersize',14,'DisplayName',[sbj.name,' ',sbj.alg])
    case 'DET'
        plot(sbj.rmse, sbj.nnz,'s','markerfacecolor',c(1,:),'markeredgecolor','k','linewidth',a,'markersize',14,'DisplayName',[sbj.name,' ',sbj.alg])
end


end

function c = getNiceColors(color_type)

dotest = false;
c1 = colormap(parula(32));
c2 = colormap(autumn(32));

if dotest
    figure('name','C1 color test');
    hold on
    for ii = 1:size(c1,1)
        plot(ii,1,'o','markerfacecolor',c1(ii,:),'markersize',12)
        text(ii-0.75,1,sprintf('%i',ii))
    end
    
    figure('name','C2 color test');
    hold on
    for ii = 1:size(c2,1)
        plot(ii,1,'o','markerfacecolor',c2(ii,:),'markersize',12)
        text(ii-0.75,1,sprintf('%i',ii))
    end
    keyboard
end

switch color_type
    case 'cold'
        c = [c1([1 3 6 9],:) ];
    case 'medium'
        c = [c1([12 16 19 23],:) ];
    case 'hot'
        c = [c2([32 25 13 5],:)];
end

end


function [sbj] = retrieve_results(fe,alg,name)
sbj.alg = alg;
sbj.name = name;

rmse = feGet(fe,'vox rmse')./feGet(fe,'b0signalimage')';
rmse = rmse(rmse~=Inf);
rmse = nanmean(rmse);
sbj.rmse = rmse;

ind = find(~isnan(fe.life.fit.weights)&fe.life.fit.weights>0);
nnzeros = length(ind);
sbj.nnz = nnzeros; 

end

