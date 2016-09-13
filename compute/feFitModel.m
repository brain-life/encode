function [fit, w, R2] = feFitModel(M,dSig,fitMethod)
% 
% feFitModel() function in LiFE but restricted to the
% 
% BBNNLS algorithm and using the Factorization model.
% M is the factorization model composed by:
%
%   M.DictSig    Dictionary
%   
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com


% Below are the old comments which are not longer correct in general
% Fit the LiFE model.
%
% Finds the weights for each fiber to best predict the directional
% diffusion signal (dSig)
%
%  fit = feFitModel(M,dSig,fitMethod)
%
% dSig:  The diffusion weighted signal measured at each
%        voxel in each direction. These are extracted from 
%        the dwi data at some white-matter coordinates.
% M:     The LiFE difusion model matrix, constructed
%        by feConnectomeBuildModel.m
%
% fitMethod: 
%  - 'bbnnls' - DEFAULT and best, faster large-scale solver.
%
% See also: feCreate.m, feConnectomeBuildModel.m, feGet.m, feSet.m
%
% Example:
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
% Notes about the LiFE model:
%
% The rows of the M matrix are nVoxels*nBvecs. We are going to predict the
% diffusion signal in each voxel for each direction.
%
% The columns of the M matrix are nFibers + nVoxels.  The diffusion signal
% for each voxel is predicted as the weighted sum of predictions from each
% fibers that passes through a voxel plus an isotropic (CSF) term.
%
% In addition to M, we typically return dSig, which is the signal measured
% at each voxel in each direction.  These are extracted from the dwi data
% and knowledge of the roiCoords.

% fit the model, by selecting the proper toolbox.

mycomputer = computer();
release = version('-release');

switch fitMethod
   case {'bbnnls'}
    [nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
    %[nTheta]  = size(M.DictSig,1);
    [nAtoms] = size(M.DictSig,2); %feGet(fe,'natoms');
    [Nvoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');
    
    tic
    fprintf('\nLiFE: Computing least-square minimization with BBNNLS...\n')
    opt = solopt;
    opt.maxit = 5000;
    opt.use_tolo = 1;
    
    switch strcat(mycomputer,'_',release)
        case {'GLNXA64_2015a'}
        out_data = bbnnls_GLNXA64(M,dSig,zeros(nFibers,1),opt);
        case {'MACI64_2014b'}
        out_data = bbnnls_MACI64(M,dSig,zeros(nFibers,1),opt);
        otherwise
        sprintf('WARNING: currently LiFE is optimized for an efficient usage of memory \n using the Sparse Tucker Decomposition aproach (Caiafa&Pestilli, 2015) \n ONLY for Linux (MatlabR2015a) and MacOS (MatlabR2014b). \n If you have a different system or version you can still \n use the old version of LiFE (memory intensive). \n\n')
        sprintf('\n Starting using old version of LiFE...\n')
        out_data = bbnnls_OLD(M.MmatrixM,dSig,zeros(nFibers,1),opt);
    end
    fprintf('BBNNLS status: %s\nReason: %s\n',out_data.status,out_data.termReason);
    w = out_data.x;
    fprintf(' ...fit process completed in %2.3fminutes\n',toc/60)
    % Save the state of the random generator so that the stochasit cfit can be recomputed.
    defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
    fit.randState = defaultStream.State;   
    
    % Save out some results 
    fit.results.R2        = [];
    fit.results.nParams   = size(M,2);
    fit.results.nMeasures = size(M,1);
    R2=[];

   otherwise
     error('Cannot fit LiFE model using method: %s.\n',fitMethod);
end

% Save output structure.
fit.weights             = w;
fit.params.fitMethod    = fitMethod;

end
