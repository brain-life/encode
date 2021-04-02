function se = feComputeEvidence_norm(rmse1, rmse2, nbins)
% Computes a series of distance metrics between two RMSE distribtions.
%
% Compute summary statistics on to characterize the lesion,:
% We compute:
% - The strength of evidence, namely the effect size of the lesion
% - The Kullback-Leibler Divergence
% - Jeffrey's Divergence
% - The Earth Mover's distance
%
% Copyright (2016), Franco Pestilli, Indiana University, frakkopesto@gmail.com.
%

% Prepare the distribution of errors and the histograms describing the errors.
se.nolesion.rmse.all  = rmse1;
se.nolesion.rmse.mean = mean(se.nolesion.rmse.all);

se.lesion.rmse.all    = rmse2;
se.lesion.rmse.mean   = mean(se.lesion.rmse.all);

% Histograms

% define the bins based on the inputs / argument
se.nbins     = nbins;
se.xrange(1) = min([se.lesion.rmse.all se.nolesion.rmse.all]);
se.xrange(2) = max([se.lesion.rmse.all se.nolesion.rmse.all]);
se.bins      = linspace(se.xrange(1),se.xrange(2),se.nbins);

% compute the histograms on the data
[se.lesion.hist,   se.lesion.xhist]   = hist(se.lesion.rmse.all, se.bins);
[se.nolesion.hist, se.nolesion.xhist] = hist(se.nolesion.rmse.all, se.bins);

% normalize the histograms
se.lesion.hist     = se.lesion.hist ./ sum(se.lesion.hist);
se.nolesion.hist   = se.nolesion.hist ./ sum(se.nolesion.hist);

%% compute the true differences between the histograms

% Kullback-Leibler Divergence
se.kl.name = sprintf('Kullback–Leibler divergence: http://en.wikipedia.org/wiki/Kullback-Leibler_divergence');
tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ (se.lesion.hist + eps) );
se.kl.mean = nansum(tmp);clear tmp
se.kl.std  = nan;

% Jeffrey's divergence
se.j.name = sprintf('Jeffrey''s divergence: http://en.wikipedia.org/wiki/Divergence_(statistics)');
tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2)  ) + ...
      se.lesion.hist .* log2( (se.lesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2) );
se.j.mean = nansum(tmp); clear tmp
se.j.std   = nan;

% Earth Mover's Distance:
% Note: This can be very slow and my require large amounts of memory for more than 1000 voxels
fprintf('[%s] Computing the Earth Mover''s distance... \n',mfilename)
se.em.name = sprintf('Earth Mover''s distance: http://en.wikipedia.org/wiki/Earth_mover''s_distance');

try % Using Rubinov c-code fastest
    pairwiseDist = zeros(size(se.lesion.xhist,2));
    for i=1:size(se.nolesion.xhist,2)
        for j=1:size(se.lesion.xhist,2)
            pairwiseDist(i,j) = abs(se.lesion.xhist(i)-se.nolesion.xhist(j));
        end
    end
    tmp_em = emd_mex(se.nolesion.hist,se.lesion.hist,pairwiseDist);
catch % else
    fprintf('[%s] Cannot find compiled c-code for Earth Movers Distance.\nUsing the slower and less reliable MatLab implementation.',mfilename)
    [~,tmp_em] = emd(se.nolesion.xhist',se.lesion.xhist',se.nolesion.hist',se.lesion.hist',@gdf);
end
se.em.mean = tmp_em;
se.em.std  = nan;
clear tmp_emp

% Strenght of evidence (effect size)
fprintf('[%s] Computing the Strength of Evidence... ', mfilename)
se.s.name = sprintf('strength of evidence, d-prime: http://en.wikipedia.org/wiki/Effect_size');
se.s.nboots = 5000; 
se.s.nmontecarlo = 5;
se.s.nbins = 200;
sizeunlesioned    = length(se.nolesion.rmse.all);
nullDistributionW = nan(se.s.nboots,se.s.nmontecarlo);
nullDistributionWO = nan(se.s.nboots,se.s.nmontecarlo);
min_x = floor(mean([se.nolesion.rmse.all]) - mean([se.nolesion.rmse.all])*.05);
max_x = ceil( mean([se.lesion.rmse.all])   + mean([se.lesion.rmse.all])*.05);

for inm = 1:se.s.nmontecarlo
    fprintf('.')
    for ibt = 1:se.s.nboots
        nullDistributionW(ibt, inm) = mean(randsample(se.nolesion.rmse.all, sizeunlesioned, true));      
        nullDistributionWO(ibt, inm) = mean(randsample(se.lesion.rmse.all, sizeunlesioned, true));
    end
    
    % build the combined null distribution
    ndist = [ se.lesion.rmse.all, se.nolesion.rmse.all ];
    
    % build the true difference to compare after null
    empDiff = se.lesion.rmse.mean - se.nolesion.rmse.mean;
    
    % for every boot strap
    for ibt = 1:se.s.nboots
        
        % pull a random sample from combined distributions
        tnullDistributionW  = mean(randsample(ndist, sizeunlesioned, true));      
        tnullDistributionWO = mean(randsample(ndist, sizeunlesioned, true));
        
        % pull the null difference for comparison
        diffDist(ibt) = tnullDistributionW - tnullDistributionWO;
        
    end
    
    % get the count of nulls that are greater than observed
    se.s.pval(inm) = sum(diffDist > empDiff) / se.s.nboots; 
    % find the code in CCA and copy it here to estimate this.
    %pval(ii) = sum(grotRr(1, ii, :) > cca.cca.hocorrs(ii)) / Nperm;    
    
    % Distribution unlesioned
    [y(:,inm),xhis] = hist(nullDistributionW(:,inm),linspace(min_x,max_x,se.s.nbins));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
    
    % Distribution lesioned
    [woy(:,inm),woxhis] = hist(nullDistributionWO(:,inm),linspace(min_x,max_x,se.s.nbins));
    woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
    
end

se.s.mean = mean(diff([mean(nullDistributionW,1); ...
                          mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
se.s.std  = std(diff([mean(nullDistributionW,1); ...
                          mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
disp(' done.')

% Extract the mean and error of the botstrapped disributions of mean errors
y_m   = mean(y,2);
ywo_m = mean(woy,2);
y_e   = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];
ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];
se.s.lesioned_e   = ywo_e;
se.s.lesioned_m   = ywo_m;
se.s.unlesioned_e = y_e;
se.s.unlesioned_m = y_m;
se.s.lesioned.xbins   = woxhis;
se.s.unlesioned.xbins = xhis;
se.s.min_x   = min_x;
se.s.max_x = max_x;

%% add the bootstrapped null / std measures for all metrics

% pull the size of a distribution
nobs1 = size(rmse1, 2);
nobs2 = size(rmse2, 2);

% build the combined null distribution
nrmse = [ rmse1, rmse2 ];

% preallocate null vector for every measure
ncor = nan(se.s.nboots, 1);
ncos = nan(se.s.nboots, 1);
nkld = nan(se.s.nboots, 1);

for ibt = 1:se.s.nboots
    
    % pull random samples from each distribution
    trmse1 = randsample(rmse1, nobs1, true); % nolesion
    trmse2 = randsample(rmse2, nobs2, true); % lesion
    
    % pull random samples from combined distributions
    nrmse1 = randsample(nrmse, nobs1, true);
    nrmse2 = randsample(nrmse, nobs2, true);
    
    % estimate each null throw on raw RMSE values
    ncor(ibt) = pdist2(nrmse1, nrmse2, 'correlation');
    ncos(ibt) = pdist2(nrmse1, nrmse2, 'cosine');
    
    % build the null histogram
    
    % define the bins based on the inputs / argument
    tse.xrange = minmax([ trmse1 trmse2 ]);
    tse.bins = linspace(tse.xrange(1), tse.xrange(2), nbins);
    
    % compute the histograms on the data
    [ tse.lesion.hist, tse.lesion.xhist ] = hist(trmse2, nbins);
    [ tse.nolesion.hist, tse.nolesion.xhist ] = hist(trmse1, nbins);
    
    % normalize the histograms
    tse.lesion.hist = tse.lesion.hist ./ sum(tse.lesion.hist);
    tse.nolesion.hist = tse.nolesion.hist ./ sum(tse.nolesion.hist);
    
    % kl divergence estimated on null
    tkld = tse.nolesion.hist .* log2( (tse.nolesion.hist) ./ (tse.lesion.hist + eps) );
    nkld(ibt) = nansum(tkld);
    
end

% dissimilarity between raw RMSE values
se.dis.name = 'Dissimilarity';
se.dis.mean = pdist2(rmse1, rmse2, 'correlation');
se.dis.pval = 1 - (sum(ncor > se.dis.mean) / se.s.nboots);

% cosine dissimilarity between raw RMSE values
se.cos.name = 'Cosine Distance';
se.cos.mean = pdist2(rmse1, rmse2, 'cosine');
se.cos.pval = 1 - (sum(ncos > se.cos.mean) / se.s.nboots);

se.nkl.name = "New KL Divergence";
se.nkl.mean = se.kl.mean;
se.nkl.pval = 1 - (sum(nkld > se.kl.mean) / se.s.nboots);

end
