function se = feComputeEMD(rmse1,rmse2)
% Computes a series of distance metrics between two RMSE distribtions.
%
% Compute summary statistics on to characterize the lesion,:
% We compute:
% - The strength of evidence, namely the effect size of the lesion
% - The Kullback-Leibler Divergence
% - Jeffrey's Divergence
% - The Eath Mover's distance
%
% Copyright (2016), Franco Pestilli, Indiana University, frakkopesto@gmail.com.

% Prepare the distribution of errors and the histograms describing the
% erros.
se.nolesion.rmse.all  = rmse1;
se.lesion.rmse.all    = rmse2;
se.nolesion.rmse.mean = mean(se.nolesion.rmse.all);
se.lesion.rmse.mean   = mean(se.lesion.rmse.all);

% Histograms
se.xrange(1) = ceil(min([se.lesion.rmse.all se.nolesion.rmse.all]));
se.xrange(2) = floor(max([se.lesion.rmse.all se.nolesion.rmse.all]));
se.nbins     = 60;
se.bins      = linspace(se.xrange(1),se.xrange(2),se.nbins);
[se.lesion.hist,   se.lesion.xhist]   = hist(se.lesion.rmse.all,  se.bins);
[se.nolesion.hist, se.nolesion.xhist] = hist(se.nolesion.rmse.all,se.bins);
se.lesion.hist     = se.lesion.hist ./   sum(se.lesion.hist);
se.nolesion.hist   = se.nolesion.hist ./ sum(se.nolesion.hist);

% % Kullback-Leibler Divergence
% se.kl.name = sprintf('Kullback–Leibler divergence: http://en.wikipedia.org/wiki/Kullback-Leibler_divergence');
% tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ (se.lesion.hist + eps) );
% se.kl.mean = nansum(tmp);clear tmp
% se.kl.std  = nan;
% 
% % Jeffrey's divergence
% se.j.name = sprintf('Jeffrey''s divergence: http://en.wikipedia.org/wiki/Divergence_(statistics)');
% tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2)  ) + ...
%       se.lesion.hist .* log2( (se.lesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2) );
% se.j.mean = nansum(tmp); clear tmp
% se.j.std   = nan;

% Earth Mover's Distance:
% Note: This can be very slow and my require large amounts of memory for more than 1000 voxels
fprintf('[%s] Computing the Earth Mover''s distance... \n',mfilename)
se.em.name = sprintf('Earth Mover''s distance: http://en.wikipedia.org/wiki/Earth_mover''s_distance');

try
%if (exist('emd_mex.m','file') == 2) % Using Rubinov c-code fastest
    pairwiseDist = zeros(size(se.lesion.xhist,2));
    for i=1:size(se.nolesion.xhist,2)
        parfor j=1:size(se.lesion.xhist,2)
            pairwiseDist(i,j) = abs(se.lesion.xhist(i)-se.nolesion.xhist(j));
        end
    end
    tmp_em = emd_mex(se.nolesion.hist,se.lesion.hist,pairwiseDist);
    disp('found Rubinov code')
catch ME %else
    fprintf('[%s] Cannot find compiled c-code for Earth Movers Distance.\nUsing the slower and less reliable MatLab implementation.',mfilename)
    [~,tmp_em] = emd(se.nolesion.xhist',se.lesion.xhist',se.nolesion.hist',se.lesion.hist',@gdf);
end
se.em.mean = tmp_em;
se.em.std  = nan;
clear tmp_emp

% % Strenght of evidence (effect size)
% fprintf('[%s] Computing the Strength of Evidence... \n',mfilename)
% se.s.name = sprintf('strength of evidence, d-prime: http://en.wikipedia.org/wiki/Effect_size');
% se.s.nboots = 5000; 
% se.s.nmontecarlo = 5;
% se.s.nbins = 200;
% sizeunlesioned    = length(se.nolesion.rmse.all);
% nullDistributionW = nan(se.s.nboots,se.s.nmontecarlo);
% nullDistributionWO = nan(se.s.nboots,se.s.nmontecarlo);
% min_x = floor(mean([se.nolesion.rmse.all]) - mean([se.nolesion.rmse.all])*.05);
% max_x = ceil( mean([se.lesion.rmse.all])   + mean([se.lesion.rmse.all])*.05);
% 
% for inm = 1:se.s.nmontecarlo
%     fprintf('.')
%     for ibt = 1:se.s.nboots
%         nullDistributionW(ibt,inm)  = mean(randsample(se.nolesion.rmse.all,   sizeunlesioned,true));      
%         nullDistributionWO(ibt,inm) = mean(randsample(se.lesion.rmse.all,sizeunlesioned,true));
%     end
%     
%     % Distribution unlesioned
%     [y(:,inm),xhis] = hist(nullDistributionW(:,inm),linspace(min_x,max_x,se.s.nbins));
%     y(:,inm) = y(:,inm)./sum(y(:,inm));
%     
%     % Distribution lesioned
%     [woy(:,inm),woxhis] = hist(nullDistributionWO(:,inm),linspace(min_x,max_x,se.s.nbins));
%     woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
% end
%     
% se.s.mean = mean(diff([mean(nullDistributionW,1); ...
%                           mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
% se.s.std  = std(diff([mean(nullDistributionW,1); ...
%                           mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
% disp(' done.')
% 
% % Extract the mean and error of the botstrapped disributions of mean errors
% y_m   = mean(y,2);
% ywo_m = mean(woy,2);
% y_e   = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];
% ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];
% se.s.lesioned_e   = ywo_e;
% se.s.lesioned_m   = ywo_m;
% se.s.unlesioned_e = y_e;
% se.s.unlesioned_m = y_m;
% se.s.lesioned.xbins   = woxhis;
% se.s.unlesioned.xbins = xhis;
% se.s.min_x   = min_x;
% se.s.max_x = max_x;

end
