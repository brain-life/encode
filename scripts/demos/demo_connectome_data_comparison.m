function [fh, fe] = demo_connectome_data_comparison()
% This demo characterizes connectomes obtained with different data sets and
% tracking methods. 
% 
% It compares two fundamental properties of a connectome density and error
% in predicting the diffusion signal. It shows how these conenctome
% properties depend from data type (spatial resolution and directional
% resolution) and tractogrpahy method (probabilistic, deterministic,
% constrained spehreical deconvolution or tensor based).
% 
% Below, we load previously computed results and show the
% relationship between the root-mean-squared error of a connectome in
% predicting the demeaned diffusion-weighted signal and the number of
% non-zero weighted connectome fascicles (connectome density).
% 
% The demo introduces uses of some fundamental operations provided by the
% toolbox and accessed via feGet.m 
%
% feGet.m is a hub function that allows computing connectome error and
% density among other operations.
%
% The plots geenrated by this demo reproduce partially results presented in
% Figure 3 of "Multidimensional encoding of brain connectomes" by Cesar F.
% Caiafa and Franco Pestilli, submitted (2016).
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa
%  (CONICET) email: frakkopesto@gmail.com and ccaiafa@gmail.com

%% (0) Check matlab, data dependencies and path settings.
if ~exist('vistaRootPath.m','file');
    disp('Vistasoft package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/vistalab/vistasoft');
end
if ~exist('feDemoDataPath.m','file');
    disp('ERROR: demo dataset either not installed or not on matlab path.')
    error('Please, download it from http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz')
end
% s = what('demo_datasets');
% if isempty(s)
%      disp('ERROR: demo dataset either not installed or not on matlab path.')
%      error('Please, download it from http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz')
% end

%% (1) Figure 3 from Multidimensional encoding of brain connectomes
%      Cesar F. Caiafa and Franco Pestilli, submitted.
%
% Below we will first load part of the data used in Figure 3 of the
% original publication (Caiafa and Pestilli, submitted)
%
% After that we will add one additional point to the plot. Using data from
% a different subject.
%
% This plots shows in a compact form two fundamental properties of a brain
% connectome:
% - the error of the conenctome in predicting the measured diffusion
%   signal, the root-mean-squared error.
% - the density of a connectome. More specifcially the number of fibers
%   supported by the measured diffusion-weighted data in the provided
%   tractography solution.
%Generate_Fig3_paper_Caiafa_Pestilli('original')
%savefig('Fig_3a_paper.fig')
openfig('Fig_3a_paper.fig')

% We brighten the symbols to use them as background.
%Generate_Fig3_paper_Caiafa_Pestilli('gray')
%savefig('Fig_3a_paper_gray.fig')
openfig('Fig_3a_paper_gray.fig')

%% (2) Read HCP3T subject connectome obtained by using Probabilistic tractography
%
% We load data not yet present on the plot.
%
disp('loading fe_structures for 105115 subject in HCP3T dataset (PROB) ...')
feFileName = fullfile(feDemoDataPath('HCP3T','sub-105115','fe_structures'), ...
            'fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat');
%feFileName = fullfile(s.path,'HCP3T','sub-105115','fe_structures', 'fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat');         
load(feFileName)

% Here we extract two measures we are interested in:
% (1) The root-mean-squared-error RMSE of the connectome in predicting 
%     the measured demeaned diffusion-weighted signal.
% (2) The number of non-zero weighted fibers. These are fibers for which
%     LiFE assiged a weight larger than zero.
%
% First we pick a data set.
sbj.alg = 'PROB';
sbj.name = 'HCP3T';

% We use the core function feGet.m to extract the RMSE and the B0 (MRI
% measureemnts without the diffusion-weighted gradient applied).
% 
% We compute the mean RMSE across the whole white matter volume.
sbj.rmse = nanmean(feGet(fe,'voxrmses0norm'));

% We find the positive weights and disregard the NaNs. THen compute the
% number of postive weights (number of fascicles with non-zero weight, alse
% referred to as conenctome density).
sbj.nnz = feGet(fe,'connectome density'); 

% Finally we add the new data point to the plot we have generted. This si
% doen by plotting connectome density on the ordinate and RMSE on the
% abscissa.
Add_new_data_point(sbj,'cold',2)
%
% Below we show additional examples of data points added to the principla
% plot in the Figure. To do so, we repeate several operations shown abouve
% (using feGet.m). But we reduce clutter by packaging the operations into a
% helper function saved at the bottom of this file that can be conveneinty
% called multiple times as we show examples of multiple data points added
% to the plot.
%

%% (3) Read data from the chosen subjects.
%
% 3.1 These results were obtained by using tensor-based deterministic
% tractography and the HCP3T data set.
%
% In practice we repeate the same operations shown about (using feGet.m).
% To reduce clutter we have packaged the operations into a helper function
% saved at the bottom of this file that can be conveneinty called multiple
% times as we show examples of multiple data points added to the plot.
%
disp('loading fe_structures for 105115 subject in HCP3T dataset (DET) ...')
feFileName = fullfile(feDemoDataPath('HCP3T','sub-105115','fe_structures'), ...
             'fe_structure_105115_STC_run01_tensor__connNUM01.mat');
%feFileName = fullfile(s.path,'HCP3T','sub-105115','fe_structures','fe_structure_105115_STC_run01_tensor__connNUM01.mat');
load(feFileName)
sbj = retrieve_results(fe,'TENSOR', 'HCP3T');

% plot new data point
Add_new_data_point(sbj,'cold',2)

% 3.2 These results were obtained by using CSD-based Probabilistic
% tractography and the STN data set.
disp('loading fe_structures for FP subject in STN dataset (PROB) ...')
 feFileName = fullfile(feDemoDataPath('STN','sub-FP','fe_structures'), ...
              'fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_PROB_lmax10_connNUM01.mat');
%feFileName = fullfile(s.path,'STN','sub-FP','fe_structures','fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_SD_PROB_lmax10_connNUM01.mat');
load(feFileName)
sbj = retrieve_results(fe,'PROB', 'STN');

% plot new data point
Add_new_data_point(sbj,'medium',2)

% 3.3 These results were obtained by using tensor-based deterministic
% tractography and the STN data set.
disp('loading fe_structures for FP subject in STN dataset (DET) ...')
 feFileName = fullfile(feDemoDataPath('STN','sub-FP','fe_structures'), ...
              'fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_tensor__connNUM01.mat');
%feFileName = fullfile(s.path,'STN','sub-FP','fe_structures', 'fe_structure_FP_96dirs_b2000_1p5iso_STC_run01_tensor__connNUM01.mat');
load(feFileName)
sbj = retrieve_results(fe,'TENSOR', 'STN');

% plot new data point
Add_new_data_point(sbj,'medium',2)

% 3.4 These results were obtained by using CSD-based probabilistic
% tractography and the HCP7T data set.
disp('loading fe_structures for 108323 subject in HCP7T dataset (PROB) ...')
 feFileName = fullfile(feDemoDataPath('HCP7T','sub-108323','fe_structures'), ...
              'fe_structure_108323_STC_run01_SD_PROB_lmax8_connNUM01.mat');
%feFileName = fullfile(s.path,'HCP7T','sub-108323','fe_structures','fe_structure_108323_STC_run01_SD_PROB_lmax8_connNUM01.mat');
load(feFileName)
sbj = retrieve_results(fe,'PROB', 'HCP7T');

% plot new data point
Add_new_data_point(sbj,'hot',2)

% 3.5 These results were obtained by using tensor-based deterministic
% tractography and the HCP7T data set.
disp('loading fe_structures for 108323 subject in HCP7T dataset (DET) ...')

 feFileName = fullfile(feDemoDataPath('HCP7T','sub-108323','fe_structures'), ...
              'fe_structure_108323_STC_run01_tensor__connNUM01.mat');
%feFileName = fullfile(s.path,'HCP7T','sub-108323','fe_structures','fe_structure_108323_STC_run01_tensor__connNUM01.mat');
load(feFileName)
sbj = retrieve_results(fe,'TENSOR', 'HCP7T');

% plot new data point
Add_new_data_point(sbj,'hot',2)


end

% Below is a series of local helper functions.
function [] = Generate_Fig3_paper_Caiafa_Pestilli(color_mode)
%
% Load data from the demo data repositroy and geenrate a plot similar to
% the one in Figure 3 of Caiafa and Pestilli under review.
%

%DataPath = feDemoDataPath('Figs_data');
DataPath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Revision_Feb2017/Results/Variability/';

HCP_subject_set = {'111312','105115','113619','110411'};
STN_subject_set = {'KK_96dirs_b2000_1p5iso','FP_96dirs_b2000_1p5iso','HT_96dirs_b2000_1p5iso','MP_96dirs_b2000_1p5iso'};
HCP7T_subject_set = {'108323','109123','111312_7T','125525','102311_Paolo_masks'};

fh = figure('name','combined scatter mean +-sem across repeats','color','w');
set(fh,'Position',[0,0,800,600]);

Nalg = 13; % We plot a few data points (13 in total, 6 Prob + 6 Stream + Tensor)

% plot HCP
Gen_plot(HCP_subject_set,'cold',DataPath,Nalg,'HCP3T90',color_mode)

% plot STN
Gen_plot(STN_subject_set,'medium',DataPath,Nalg,'STN96',color_mode)

Nalg = 9; % We plot a few data points (9 in total, 4 Prob + 4 Stream + Tensor)

% plot HCP7T
Gen_plot(HCP7T_subject_set,'hot',DataPath,Nalg,'HCP7T60',color_mode)

set(gca,'tickdir','out', 'ticklen',[0.025 0.025], ...
         'box','off','ytick',[2 15 32].*10^4, 'xtick', [0.04 0.07 0.1], ...
         'ylim',[2 32].*10^4, 'xlim', [0.04 0.1],'fontsize',20)
axis square
ylabel('Fascicles number','fontsize',20)
xlabel('Connectome error (r.m.s.)','fontsize',20)
drawnow

end

function [] = Gen_plot(subject_set,color_type,DataPath,Nalg,dataset,color_mode)
%
% Generate a scatter plot similar to Caiafa and Pestilli Figure 3
%
nnz_all = zeros(length(subject_set),Nalg,10);
nnz_mean = zeros(length(subject_set),Nalg);
nnz_std  = zeros(length(subject_set),Nalg);

alg_names = cell(1,Nalg);

if Nalg==13
    range_prob = 2:2:12;
    range_det = 3:2:13;
    prob_ix_low = [2:7];
    prob_ix_high = [15:20];
    det_ix_low = [8:13];
    det_ix_high = [21:26];
    ten_ix_low = [1];
    ten_ix_high = [14];
    lmax_order = [3,4,5,6,1,2];
else
    range_prob = 2:2:8;
    range_det = 3:2:9;
    prob_ix_low = [2:5];
    prob_ix_high = [11:14];
    det_ix_low = [6:9];
    det_ix_high = [15:18];
    ten_ix_low = [1];
    ten_ix_high = [10];
    lmax_order = [1,2,3,4];
end

n = 1;
for subject = subject_set;
    switch dataset
        case {'HCP7T60','STN96','HCP3T90'}
            DataFile = char(fullfile(DataPath,strcat('Rmse_nnz_10_connectomes_',subject,'_run01','.mat')));
        case {'HCP3T60','STN60'}
            DataFile = char(fullfile(DataPath,strcat('Rmse_nnz_10_connectomes_',subject,'_60dir*run01','.mat'))); 
    end    
    
    load(DataFile)
    
    m = 1;
    % Tensor
    for p=1:1
        rmse_all(n,m,:) = Result_alg(p).rmse;
        rmse_mean(n,m)  = nanmean(Result_alg(p).rmse);
        rmse_std(n,m)   = nanstd(Result_alg(p).rmse)./sqrt(length(Result_alg(p).rmse));
        
        nnz_all(n,m,:) = Result_alg(p).nnz;
        nnz_mean(n,m)  = nanmean(Result_alg(p).nnz);
        nnz_std(n,m)   = nanstd(Result_alg(p).nnz)./sqrt(length(Result_alg(p).nnz));
        
        alg_names{m} = char(alg_info(p).description);
        m = m +1;
    end

    % Prob
    for p = range_prob  
        rmse_all(n,m,:) = Result_alg(p).rmse;
        rmse_mean(n,m) = nanmean(Result_alg(p).rmse);
        rmse_std(n,m)  = nanstd(Result_alg(p).rmse)./sqrt(length(Result_alg(p).rmse));       
        
        nnz_all(n,m,:) = Result_alg(p).nnz;
        nnz_mean(n,m)  = mean(Result_alg(p).nnz);
        nnz_std(n,m)   = std(Result_alg(p).nnz)./sqrt(length(Result_alg(p).nnz));
        
        alg_names{m} = char(alg_info(p).description);
        m = m +1;
    end

    % Det
    for p = range_det           
        rmse_all(n,m,:) = Result_alg(p).rmse;
        rmse_mean(n,m) = nanmean(Result_alg(p).rmse);
        rmse_std(n,m)  = nanstd(Result_alg(p).rmse)./sqrt(length(Result_alg(p).rmse));
        
        nnz_all(n,m,:) = Result_alg(p).nnz;
        nnz_mean(n,m) = nanmean(Result_alg(p).nnz);
        nnz_std(n,m)  = nanstd(Result_alg(p).nnz)./sqrt(length(Result_alg(p).nnz)); 
        
        alg_names{m} = char(alg_info(p).description);
        m = m +1;
    end

    n = n + 1;
end

switch color_mode
    case 'original'
        c = getNiceColors(color_type);
    case 'gray'
        c = repmat([.9,.9,.9], [length(subject_set),1]);
end


for is  = 1:size(nnz_all,1)    
    tmp_rmse = squeeze(rmse_all(is,:,:));
    tmp_rmse(isinf(tmp_rmse)) = nan;
    
    tmp_nnz = squeeze(nnz_all(is,:,:));
    tmp_nnz(isinf(tmp_nnz)) = nan;
    
    % mu and sem RMSE
    rmse_mu(is,:) = squeeze(nanmean(tmp_rmse,2));
    rmse_ci(is,:) = [rmse_mu(is,:), rmse_mu(is,:)] + 5*([-nanstd(tmp_rmse,[],2),;nanstd(tmp_rmse,[],2)]' ./sqrt(size(tmp_rmse,2)));
    
    % mu and sem NNZ
    nnz_mu(is,:) = squeeze(nanmean(tmp_nnz,2));
    nnz_ci(is,:) = [nnz_mu(is,:), nnz_mu(is,:)] + 5*([-nanstd(tmp_nnz,[],2);nanstd(tmp_nnz,[],2)]' ./sqrt(size(tmp_rmse,2)));
end

% scatter plot with confidence intervals first all in gray
a = 0.5;

for ii = 1:length(subject_set) % subjects
   hold on
   % PROB
   for iii = lmax_order
       plot(rmse_mean(ii,prob_ix_low(iii)), nnz_mean(ii,prob_ix_low(iii)),'o','markerfacecolor',c(ii,:),'markeredgecolor',[.5,.5,.5],'linewidth',0.5,'markersize',14)
       plot([rmse_ci(ii,prob_ix_low(iii)); rmse_ci(ii,prob_ix_high(iii))], [nnz_mu(ii,prob_ix_low(iii)); nnz_mu(ii,prob_ix_low(iii))],'-','color',[a a a],'linewidth',2)
       plot([rmse_mu(ii,prob_ix_low(iii)); rmse_mu(ii,prob_ix_low(iii))], [nnz_ci(ii,[prob_ix_low(iii)]);  nnz_ci(ii,prob_ix_high(iii))],'-','color',[a a a],'linewidth',2)   
   end
   
   % DET
   for iii = lmax_order
       plot(rmse_mean(ii,det_ix_low(iii)), nnz_mean(ii,det_ix_low(iii)),'s','markerfacecolor',c(ii,:),'markeredgecolor',[.5,.5,.5],'linewidth',0.5,'markersize',14)
       plot([rmse_ci(ii,det_ix_low(iii)); rmse_ci(ii,[det_ix_high(iii)])], [nnz_mu(ii,det_ix_low(iii)); nnz_mu(ii,det_ix_low(iii))],'-','color',[a a a],'linewidth',2)
       plot([rmse_mu(ii,det_ix_low(iii)); rmse_mu(ii,det_ix_low(iii));], [nnz_ci(ii,det_ix_low(iii)); nnz_ci(ii,[det_ix_high(iii)])],'-','color',[a a a],'linewidth',2)
   end
   
   % TENSOR
   plot(rmse_mean(ii,ten_ix_low), nnz_mean(ii,ten_ix_low),'d','markerfacecolor',c(ii,:),'markeredgecolor',[.5,.5,.5],'linewidth',0.5,'markersize',14)
   plot([rmse_ci(ii,ten_ix_low); rmse_ci(ii,ten_ix_high)], [nnz_mu(ii,ten_ix_low); nnz_mu(ii,ten_ix_low)],'-','color',[a a a],'linewidth',2)
   plot([rmse_mu(ii,ten_ix_low); rmse_mu(ii,ten_ix_low)], [nnz_ci(ii,ten_ix_low); nnz_ci(ii,ten_ix_high)],'-','color',[a a a],'linewidth',2)
   
end

end


function c = getNiceColors(color_type)
%
% Load look-up-table for plot colors.
% 
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
        %c = [c2([32 25 13 5],:)];
        c = [c2([32 27 19 12 2],:)];
end

end

function [sbj] = retrieve_results(fe,alg,name)
%
% Extracts results from a precomputed FE strcuture.
% These results compare
%
sbj.alg = alg;
sbj.name = name;

% We use the core function feGet.m to extract the RMSE and the B0 (MRI
% measureemnts without the diffusion-weighted gradient applied).
% 
% We compute the mean RMSE across the whole white matter volume.
sbj.rmse = nanmean(feGet(fe,'voxrmses0norm'));

% We find the positive weights and disregard the NaNs. THen compute the
% number of postive weights (number of fascicles with non-zero weight, alse
% referred to as conenctome density).
sbj.nnz = feGet(fe,'connectome density'); 

end


function [] = Add_new_data_point(sbj,color_type,order)
%
% This function adds a new data point precomputed into the scater plot that
% compares connectome prediction error and resolution.
%
c = getNiceColors(color_type);

%% scatter plot
a = 0.5;

switch sbj.alg
    case 'PROB'
        plot(sbj.rmse, sbj.nnz,'o','markerfacecolor',c(order,:),'markeredgecolor','k','linewidth',a,'markersize',14,'DisplayName',[sbj.name,' ',sbj.alg])
    case 'DET'
        plot(sbj.rmse, sbj.nnz,'s','markerfacecolor',c(order,:),'markeredgecolor','k','linewidth',a,'markersize',14,'DisplayName',[sbj.name,' ',sbj.alg])
    case 'TENSOR'
        plot(sbj.rmse, sbj.nnz,'d','markerfacecolor',c(order,:),'markeredgecolor','k','linewidth',a,'markersize',14,'DisplayName',[sbj.name,' ',sbj.alg])
end

drawnow

end




