function [ fns, pval, fdat ] = feFiberNetworkStats( data, nreps, prop, kcore )
%% given a connection matrix, returns mutliple netowrk evaluation metrics
%
%   needs BCT / NBS (?)

% Try to add:
% - gtom
% - edge_nei_overlap_bu
% - rentian_scaling (need XYZ of coords)
% - assortativity
% - local_assortativity_wu
% - participation_coef
% - gateway_coef
% - diveristy_coef

% - motifs?
% - lattice null? D is distance matrix
% - write fxn to generate plot of network parameters at each iteration of
% randomization to see if 10 is enough (Olaf's suggestion)

%% fiber data - connection matrix

% has raw metric, proportional binarized at min threshold, and normalized thresheld metric
% normalize on thresheld or raw data?

display('Thresholding Data...');

fdat.prp = prop;
fdat.kcore = kcore;
fdat.raw = data;
fdat.thr = threshold_proportional(fdat.raw, fdat.prp);
fdat.nrm = weight_conversion(fdat.thr, 'normalize');
fdat.len = weight_conversion(fdat.thr, 'lengths');
fdat.bin = weight_conversion(fdat.thr, 'binarize');
[ fdat.dist, fdat.edge ] = distance_wei(fdat.len); 
%[ fdat.tree, fdat.plus ] = backbone_wu(fdat.nrm, 5);

%% compute neighborhood metric

display('Computing Louvain Community Structure...');

for ii = 1:nreps
    % estimate community estimates on input data
    [ ci(:, ii), q(ii) ] = community_louvain(fdat.nrm, 1, [], 'modularity');
end;

% take the highest q-vaule (community structure statistic)
% - find the highest estimated number of neighborhoods
% - sort by neighborhoods
[ ~, indx_max ] = max(q);

% group nodes
ci_max = ci(:, indx_max);

% sort matrix by nodes
[ ~, bb ] = sort(ci_max);

% overlap of nieghborhood assignments
agree = agreement(ci);

% normalize for resampling N
agree = agree ./ nreps;

fdat.agree = agree(bb, bb);
fdat.nsrta = fdat.nrm(bb, bb);

% build confidend
tau = 0.5;
fns.nrm.cnsns = consensus_und(agree, tau, nreps);

%% fiber network statistics

% local
% - strenth: *.str
% - efficientcy: *.lcEff
% - betweenness centrality: *.bcv

% global
% - characteristic path length: *.chpl
% - mean clustering coefficient: mean(*.ccoef) 
% - small worldness - requires NBS / null model comparisons: 

% Daianu 2013
% - k-core

%% stats calculation across network thresholds

display('Computing Measures on Data Matrix...');

% stats on normalized input data matrix

% node measures
fns.nrm.deg = degrees_und(fdat.nrm)';
fns.nrm.str = strengths_und(fdat.nrm)';
fns.nrm.btw = betweenness_wei(fdat.dist);
fns.nrm.lcEff = efficiency_wei(fdat.nrm, 1);
fns.nrm.ccoef = clustering_coef_wu(fdat.nrm);
[ fdat.btcm, fns.nrm.bcv ] = edge_betweenness_wei(fdat.dist);
fns.nrm.eigv = eigenvector_centrality_und(fdat.nrm);

% community structure - local / global measures
[ fns.nrm.asgn, fns.nrm.stat ] = community_louvain(fdat.nrm, 1, [], 'modularity');
[ fns.nrm.struc, fns.nrm.mod ] = modularity_und(fdat.nrm, 1);
fns.nrm.pcoef = participation_coef(fdat.nrm, fns.nrm.asgn);
fns.nrm.modz = module_degree_zscore(fdat.nrm, fns.nrm.asgn);

% global measures
fns.nrm.mcoef = mean(fns.nrm.ccoef);
[ fns.nrm.dens, fns.nrm.dvrt, fns.nrm.dedg ] = density_und(fdat.nrm);
fns.nrm.trans = transitivity_wu(fdat.nrm);
fns.nrm.glEff = efficiency_wei(fdat.nrm, 0);
[ fns.nrm.chpl, ~ ] = charpath(fdat.nrm);
fns.nrm.rcc = rich_club_wu(fdat.nrm);

% stats on binarized matrix / k-core
[ fdat.kcore, fns.bin.kn, fns.bin.plord, fns.bin.pllvl ] = kcore_bu(fdat.bin, fdat.kcore);
[ fns.bin.core, fns.bin.ckn ] = kcoreness_centrality_bu(fdat.bin);

% reorient kcore stats
fns.bin.core = fns.bin.core';
fns.bin.ckn = fns.bin.ckn';

%% run null comparison test to repeated generation of random networks

display(['Running ', num2str(nreps), ' Replications on Randomized Input Matrix for Null Test Comparison...']);
for ii = 1:nreps
    
    % randomize matrix - create normalization for stats
    rmat = randmio_und_connected(fdat.nrm, 10);
    rlen = weight_conversion(rmat, 'lengths');
    rbin = weight_conversion(rmat, 'binarize');
    rdst = distance_wei(rlen);
    
    % pull and keep all network values for null comparison
    
    % node measures
    rep.deg{ii} = degrees_und(rmat)';
    rep.str{ii} = strengths_und(rmat)';
    rep.btw{ii} = betweenness_wei(rdst);
    rep.lcEff{ii} = efficiency_wei(rmat, 1);
    rep.ccoef{ii} = clustering_coef_wu(rmat);
    [ ~, rep.bcv{ii} ] = edge_betweenness_wei(rdst);
    rep.eigv{ii} = eigenvector_centrality_und(rmat);
    
    % community structure - local / global measures
    [ rep.asgn{ii}, rep.stat{ii} ] = community_louvain(rmat, 1, [], 'modularity');
    [ rep.struc{ii}, rep.mod{ii} ] = modularity_und(rmat, 1);
    rep.pcoef{ii} = participation_coef(rmat, rep.asgn{ii});
    rep.modz{ii} = module_degree_zscore(rmat, rep.asgn{ii});
    
    % global measures
    rep.mcoef{ii} = mean(rep.ccoef{ii});
    [ rep.dens{ii}, rep.dvrt{ii}, rep.dedg{ii} ] = density_und(rmat);
    rep.trans{ii} = transitivity_wu(rmat);
    rep.glEff{ii} = efficiency_wei(rmat, 0);
    [ rep.chpl{ii}, ~ ] = charpath(rmat);
    rep.rcc{ii} = rich_club_wu(rmat);

%     % stats on binarized matrix / k-core
%     [ ~, rep.kn{ii}, rep.plord{ii}, rep.pllvl{ii} ] = kcore_bu(rbin, fdat.kcore);
%     [ rep.core{ii}, rep.ckn{ii} ] = kcoreness_centrality_bu(rbin);
%     rep.core{ii} = rep.core{ii}';
%     rep.ckn{ii} = rep.ckn{ii}';
    
end

display('Calculating p-values...');

% pull counts of random network estimate exceeding data estimate
for ii = 1:nreps
    
    % node measures
    pv.deg{ii} = rep.deg{ii} >= fns.nrm.deg;
    pv.str{ii} = rep.str{ii} >= fns.nrm.str;
    pv.btw{ii} = rep.btw{ii} >= fns.nrm.btw;
    pv.lcEff{ii} = rep.lcEff{ii} >= fns.nrm.lcEff;
    pv.ccoef{ii} = rep.ccoef{ii} >= fns.nrm.ccoef;
    pv.bcv{ii} = rep.bcv{ii} >= fns.nrm.bcv;
    pv.eigv{ii} = rep.eigv{ii} >= fns.nrm.eigv;
    
    % community structure - local / global measures
    pv.asgn{ii} = rep.asgn{ii} >= fns.nrm.asgn;
    pv.stat{ii} = rep.stat{ii} >= fns.nrm.stat;
    pv.struc{ii} = rep.struc{ii} >= fns.nrm.struc;
    pv.mod{ii} = rep.mod{ii} >= fns.nrm.mod;
    pv.pcoef{ii} = rep.pcoef{ii} >= fns.nrm.pcoef;
    pv.modz{ii} = rep.modz{ii} >= fns.nrm.modz;

    % global measures
    pv.mcoef{ii} = rep.mcoef{ii} >= fns.nrm.mcoef;
    pv.dens{ii} = rep.dens{ii} >= fns.nrm.dens;
    pv.dvrt{ii} = rep.dvrt{ii} >= fns.nrm.dvrt;
    pv.dedg{ii} = rep.dedg{ii} >= fns.nrm.dedg;
    pv.trans{ii} = rep.trans{ii} >= fns.nrm.trans;
    pv.glEff{ii} = rep.glEff{ii} >= fns.nrm.glEff;
    pv.chpl{ii} = rep.chpl{ii} >= fns.nrm.chpl;
    pv.rcc{ii} = rep.rcc{ii} >= fns.nrm.rcc;
    
%     % stats on binarized matrix / k-core
%     pv.kn{ii} = rep.kn{ii} >= fns.bin.kn;
%     pv.plord{ii} = rep.plord{ii} >= fns.bin.plord;
%     pv.pllvl{ii} = rep.pllvl{ii} >= fns.bin.pllvl;
%     pv.core{ii} = rep.core{ii} >= fns.bin.core;
%     pv.ckn{ii} = rep.ckn{ii} >= fns.bin.ckn;
    
end

% sum counts for p-value

% node measures
pval.deg = sum(cell2mat(pv.deg), 2) / nreps;
pval.str = sum(cell2mat(pv.str), 2) / nreps;
pval.btw = sum(cell2mat(pv.btw), 2) / nreps;
pval.lcEff = sum(cell2mat(pv.lcEff), 2) / nreps;
pval.ccoef = sum(cell2mat(pv.ccoef), 2) / nreps;
pval.bcv = sum(cell2mat(pv.bcv), 2) / nreps;
pval.eigv = sum(cell2mat(pv.eigv), 2) / nreps;

% community structure - local / global measures
pval.asgn = sum(cell2mat(pv.asgn), 2) / nreps;
pval.stat = sum(cell2mat(pv.stat), 2) / nreps;
pval.struc = sum(cell2mat(pv.struc), 2) / nreps;
pval.mod = sum(cell2mat(pv.mod), 2) / nreps;
pval.pcoef = sum(cell2mat(pv.pcoef), 2) / nreps;
pval.modz = sum(cell2mat(pv.modz), 2) / nreps;

% global measures
pval.mcoef = sum(cell2mat(pv.mcoef), 2) / nreps;
pval.dens = sum(cell2mat(pv.dens), 2) / nreps;
pval.dvrt = sum(cell2mat(pv.dvrt), 2) / nreps;
pval.dedg = sum(cell2mat(pv.dedg), 2) / nreps; 
pval.trans = sum(cell2mat(pv.trans), 2) / nreps;
pval.glEff = sum(cell2mat(pv.glEff), 2) / nreps;
pval.chpl = sum(cell2mat(pv.chpl), 2) / nreps;
pval.rcc = sum(cell2mat(pv.rcc), 2) / nreps;

% % stats on binarized matrix / k-core
% pval.kn = rep.kn{ii} >= fns.bin.kn;
% pval.plord = rep.plord{ii} >= fns.bin.plord;
% pval.pllvl = rep.pllvl{ii} >= fns.bin.pllvl;
% pval.core = rep.core{ii} >= fns.bin.core;
% pval.ckn = rep.ckn{ii} >= fns.bin.ckn;

% find pval:
% sum(rep.stat(:) >= fns.*) / nreps;

% small-worldness calculation
sw1 = mean(cellfun(@mean, rep.mcoef));
sw2 = mean(cellfun(@mean, rep.chpl));

fns.nrm.smwrld = (fns.nrm.mcoef / sw1) / (fns.nrm.chpl / sw2);

% %% create simplest output object of network stats for now
% 
% % node metrics
% out.nodes.strength.data = fns.nrm.str;
% out.nodes.strength.pval = pval.str;
% 
% out.nodes.efficiency.data = fns.nrm.lcEff;
% out.nodes.efficiency.pval = pval.lcEff;
% 
% out.nodes.btw_center.data = fns.nrm.bcv;
% out.nodes.btw_center.pval = pval.bcv;
% 
% % global metrics
% out.global.char_pathl.data = fns.nrm.chpl;
% out.global.char_pathl.pval = pval.chpl;
% 
% out.global.mean_ccoef.data = fns.nrm.mcoef;
% out.global.mean_ccoef.pval = pval.mcoef;
% 
% out.global.smwrld_unc.data = out.global.mean_ccoef.data / out.global.char_pathl.data;
% out.global.smallworld.data = fns.nrm.smwrld;
% 
% % k-core calculation
% out.kcore.matrix = fdat.kcore;
% out.kcore.nnodes = fns.bin.kn;
% 
end

