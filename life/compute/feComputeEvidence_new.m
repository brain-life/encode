function se = feComputeEvidence_new(rmse1, rmse2, nbins, nmc, nboots)
% Computes a series of distance metrics between two RMSE distribtions.
%
% Compute summary statistics on to characterize the lesion,:
% We compute:
% - The strength of evidence, namely the effect size of the lesion
% - The Kullback-Leibler Divergence
% - Jeffrey's Divergence
% - The Earth Mover's distance
%
% Copyright (2020), Brent McPherson & Franco Pestilli, Indiana University
%

% fill in default parameters

if(~exist('nbins', 'var') || isempty(nbins))
    nbins = 128;
end

if(~exist('nmc', 'var') || isempty(nmc))
    nmc = 5;
end

if(~exist('nboots', 'var') || isempty(nboots))
    nboots = 250; % this almost certainly has to be higher, low for test
end

% pull the size of a distribution
nobs1 = size(rmse1, 2); % nolesion
nobs2 = size(rmse2, 2); % lesion

% build the combined null distribution
nrmse = [ rmse1, rmse2 ];

% the min/max observed values for SOE / d-prime histogram computation
minx = floor(mean(rmse1) - mean(rmse1)*.05);
maxx = ceil(mean(rmse2) + mean(rmse2)*.05);

% prepare the distribution of errors and the histograms describing the errors.
se.nolesion.rmse.all = rmse1;
se.nolesion.rmse.mean = mean(se.nolesion.rmse.all);

se.lesion.rmse.all = rmse2;
se.lesion.rmse.mean = mean(se.lesion.rmse.all);

% compute histograms

% define the bins based on the inputs / argument
se.nbins = nbins;
se.xrange = minmax([ se.lesion.rmse.all se.nolesion.rmse.all ]);
se.bins = linspace(se.xrange(1), se.xrange(2), se.nbins);

% store params for reference
se.params.nmontecarlo = nmc;
se.params.nboots = nboots;

% compute the histograms on the data
% the se.bins center have to be wrapped with inf to actually match hist
% https://www.mathworks.com/matlabcentral/answers/377930-histcounts-error-in-place-of-histc
[ se.lesion.hist, se.lesion.xhist ] = histcounts(se.lesion.rmse.all, [ se.bins inf ], 'Normalization', 'probability');
[ se.nolesion.hist, se.nolesion.xhist ] = histcounts(se.nolesion.rmse.all, [ se.bins inf ], 'Normalization', 'probability');

% drop inf from end to keep bin centers the same size, b/c this is "fixed"
se.lesion.xhist = se.lesion.xhist(1:end-1);
se.nolesion.xhist = se.nolesion.xhist(1:end-1);

%% compute the true differences between the rmse / histograms

% dissimilarity between raw RMSE values
se.dis.name = 'Dissimilarity: https://nikokriegeskorte.org/category/representational-similarity-analysis/';
se.dis.mean = pdist2(rmse1, rmse2, 'correlation');

% cosine dissimilarity between raw RMSE values
se.cos.name = 'Cosine Similarity: https://en.wikipedia.org/wiki/Cosine_similarity';
se.cos.mean = pdist2(rmse1, rmse2, 'cosine');

% Kullback-Leibler Divergence
se.kl.name = sprintf('Kullback–Leibler divergence: http://en.wikipedia.org/wiki/Kullback-Leibler_divergence');
tkl = se.nolesion.hist .* log2( (se.nolesion.hist) ./ (se.lesion.hist + eps) );
se.kl.mean = nansum(tkl);
clear tkl

% Jeffrey's divergence
se.j.name = sprintf('Jeffrey''s divergence: http://en.wikipedia.org/wiki/Divergence_(statistics)');
tjd = se.nolesion.hist .* log2( (se.nolesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2)  ) + ...
      se.lesion.hist .* log2( (se.lesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2) );
se.j.mean = nansum(tjd); clear tjd

% Earth Mover's Distance:
% Note: This can be very slow and my require large amounts of memory for more than 1000 voxels
se.em.name = sprintf('Earth Mover''s distance: http://en.wikipedia.org/wiki/Earth_mover''s_distance');

try % Using Rubinov c-code fastest
    
    % estimate the distance between the bin centers
    eDist = pdist2(se.nolesion.xhist', se.lesion.xhist', @(xi, xj) abs(xi-xj));
    
    % compute the emd
    tmp_em = emd_mex(se.nolesion.hist, se.lesion.hist, eDist);

catch % else

    % compute EMD using slower MATLAB fxn
    [ ~, tmp_em ] = emd(se.nolesion.xhist',se.lesion.xhist',se.nolesion.hist',se.lesion.hist',@gdf);
    
end
se.em.mean = tmp_em;
clear tmp_emp

%% add the bootstrapped null / se measures for all metrics

% define null histograms for SOE estimate
mhist1 = nan(nbins, nmc);
mhist2 = nan(nbins, nmc);

% preallocate test vectors for every measure (SE)
tcor = nan(nboots, nmc);
tcos = nan(nboots, nmc);
tkld = nan(nboots, nmc);
tjfd = nan(nboots, nmc);
temd = nan(nboots, nmc);
tdp1 = nan(nboots, nmc);
tdp2 = nan(nboots, nmc);

% preallocate null vectors for every measure (pval)
ncor = nan(nboots, nmc);
ncos = nan(nboots, nmc);
nkld = nan(nboots, nmc);
njfd = nan(nboots, nmc);
nemd = nan(nboots, nmc);
ndpd = nan(nboots, nmc);

for mc = 1:nmc
    fprintf('.')
    for boot = 1:nboots
        
        % pull random samples from each distribution for SE estimates
        trmse1 = randsample(rmse1, nobs1, true); % nolesion
        trmse2 = randsample(rmse2, nobs2, true); % lesion
        
        % pull random samples from combined distributions for p-value estimates
        nrmse1 = randsample(nrmse, nobs1, true); % nolesion
        nrmse2 = randsample(nrmse, nobs2, true); % lesion
        
        % build the boostrapped histograms for test and null

        % compute the test histograms w/ the estimated bins
        [ t1hist, t1xhist ] = histcounts(trmse1, [ se.bins inf ], 'Normalization', 'probability');
        [ t2hist, t2xhist ] = histcounts(trmse2, [ se.bins inf ], 'Normalization', 'probability');
        
        % fix the "fixed" bins
        t1xhist = t1xhist(1:end-1);
        t2xhist = t2xhist(1:end-1);
        
        % build the the combined null histogram
        
        % compute the null histograms w/ the estimated bins
        [ n1hist, n1xhist ] = histcounts(nrmse1, [ se.bins inf ], 'Normalization', 'probability');
        [ n2hist, n2xhist ] = histcounts(nrmse2, [ se.bins inf ], 'Normalization', 'probability');

        % fix the "fixed" bins
        n1xhist = n1xhist(1:end-1);
        n2xhist = n2xhist(1:end-1);

        % estimate each SE / null throw on raw RMSE values

        % SOE / d-prime test averages
        tdp1(boot, mc) = mean(trmse1);
        tdp2(boot, mc) = mean(trmse2);
        
        % SOE / d-prime null difference
        ndp1 = mean(nrmse1);
        ndp2 = mean(nrmse2);
        ndpd(boot, mc) = ndp1 - ndp2;
                
        % estimate the test / null dissimilarity between rmse       
        tcor(boot, mc) = pdist2(trmse1, trmse2, 'correlation');
        ncor(boot, mc) = pdist2(nrmse1, nrmse2, 'correlation');
        
        % estimate the test / null cosine distance between rmse
        tcos(boot, mc) = pdist2(trmse1, trmse2, 'cosine');
        ncos(boot, mc) = pdist2(nrmse1, nrmse2, 'cosine');
        
        % build the SE / null estimates on the resampled / null histograms
        
        % kl divergence
        tkld(boot, mc) = nansum(t1hist .* log2( (t1hist) ./ (t2hist + eps) ));
        nkld(boot, mc) = nansum(n1hist .* log2( (n1hist) ./ (n2hist + eps) ));
        
        % jeffery's divergence
        tjfd(boot, mc) = nansum(t1hist .* log2( (t1hist) ./ ((t2hist + t2hist + eps)./2)  ) + ...
            t2hist .* log2( (t2hist) ./ ((t2hist + t2hist + eps)./2) ));
        njfd(boot, mc) = nansum(n1hist .* log2( (n1hist) ./ ((n2hist + n2hist + eps)./2)  ) + ...
            n2hist .* log2( (n2hist) ./ ((n2hist + n2hist + eps)./2) ));
        
        % Earth Mover's Distance; cased to failover to matlab fxb if MEX fxn fails
        try
            
            % estimate EMD on resampled / null
            temd(boot, mc) = emd_mex(t1hist, t2hist, eDist);
            nemd(boot, mc) = emd_mex(n1hist, n2hist, eDist);
            
        % pass to slower full matlab fxn if MEX fails
        catch
            
            % slower MATLAB fxn
            [ ~, temd(boot, mc) ] = emd(t1xhist', t2xhist', t1hist', t2hist', @gdf);
            [ ~, nemd(boot, mc) ] = emd(n1xhist', n2xhist', n1hist', n2hist', @gdf);
            
        end
        
    end
    
    % pull the estimated bins for SOE histogram
    mbins = linspace(minx, maxx, nbins);
    
    % compute and normalize histogram of unlesioned rmse
    [ mhist1(:, mc), mxhist1 ] = histcounts(tdp1(:, mc), [ mbins inf ], 'Normalization', 'probability');
    [ mhist2(:, mc), mxhist2 ] = histcounts(tdp2(:, mc), [ mbins inf ], 'Normalization', 'probability');
    
    % fix the "fixed" bins
    mxhist1 = mxhist1(1:end-1);
    mxhist2 = mxhist2(1:end-1);
    
end
disp(' done.')

%% store the permutation resluts in output struct

% permuted dissimilarity measure
se.dis.se = mean(std(tcor, [], 1));
se.dis.pval = min(sum(ncor > se.dis.mean) / nboots);

% permuted cosine distance measure
se.cos.se = mean(std(tcos, [], 1));
se.cos.pval = min(sum(ncos > se.cos.mean) / nboots);

% permuted kl values
se.kl.se = mean(std(tkld, [], 1));
se.kl.pval = min(sum(nkld > se.kl.mean) / nboots);

% permuted jd values
se.j.se = mean(std(tjfd, [], 1));
se.j.pval = min(sum(njfd > se.j.mean) / nboots);

% permuted emd values
se.em.se = mean(std(temd, [], 1));
se.em.pval = min(sum(nemd > se.em.mean) / nboots);

% pull the tested SOE mean / std across MCs
smn = diff([ mean(tdp1, 1); mean(tdp2, 1) ]);
sse = sqrt(sum([ std(tdp1, [], 1); std(tdp2, [], 1) ].^2, 1));

% collapse SOE across MC runs
se.s.mean = mean(smn./sqrt(sse).^2);
se.s.std = std(smn./sqrt(sse).^2);
se.s.pval = min(sum(ndpd > smn) ./ nboots);

% build mean histogram of SOE distributions
mn_mhist1 = mean(mhist1, 2);
mn_mhist2 = mean(mhist2, 2);

% add 2 std around mean histogram estimates (confidence itervals)
ci_mhist1 = [mn_mhist1, mn_mhist1] + 2*[-std(mhist1, [], 2), std(mhist1, [], 2)];
ci_mhist2 = [mn_mhist2, mn_mhist2] + 2*[-std(mhist2, [], 2), std(mhist2, [], 2)];

% store lesioned / unlesioned histograms
se.s.min_x = minx;
se.s.max_x = maxx;
se.s.unlesioned_mn = mn_mhist1;
se.s.unlesioned_ci = ci_mhist1;
se.s.unlesioned_xbins = mxhist1';
se.s.lesioned_mn = mn_mhist2;
se.s.lesioned_ci = ci_mhist2;
se.s.lesioned_xbins = mxhist2';

end
