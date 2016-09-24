function [fit, w, R2] = feFitModel(varargin)
M = varargin{1};
dSig = varargin{2};
fitMethod = varargin{3};
Niter = varargin{4};
preconditioner = varargin{5};


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


% Check for mexfiles and generate them if necessary
% Mtransp_times_b function
checkMexCompiled('-largeArrayDims', '-output', 'Mtransp_times_b', '-DNDEBUG','Mtransp_times_b.c', 'Mtransp_times_b_sub.c')
% M_times_w function
checkMexCompiled('-largeArrayDims', '-output', 'M_times_w', '-DNDEBUG', 'M_times_w.c', 'M_times_w_sub.c')
% compute_diag function
checkMexCompiled('-largeArrayDims', '-output', 'compute_diag', '-DNDEBUG', 'compute_diag.c', 'compute_diag_sub.c')

if nargin <6 % no initial w0 is provided
    [nFibers] = size(M.Phi,3);
    w0 = zeros(nFibers,1);
else
    w0 = varargin{6};
end


if strcmp(preconditioner,'preconditioner') 
    [nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
    h = compute_diag(M.Phi.subs(:,1), M.Phi.subs(:,3), M.Phi.vals, M.DictSig,nFibers);
    vals = M.Phi.vals./h(M.Phi.subs(:,3));
    M.Phi = sptensor(M.Phi.subs,vals,size(M.Phi));
end


switch fitMethod
   case {'bbnnls'}
    [nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
    %[nTheta]  = size(M.DictSig,1);
    [nAtoms] = size(M.DictSig,2); %feGet(fe,'natoms');
    [Nvoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');
    
    tic
    fprintf('\nLiFE: Computing least-square minimization with BBNNLS...\n')
    opt = solopt;
    opt.maxit = Niter;
    opt.use_tolo = 1;
    opt.tolg = 1e-5;
    opt.verbose = 1;
    
    out_data = bbnnls(M,dSig,w0,opt);
    
    if strcmp(preconditioner,'preconditioner')
        out_data.x = out_data.x./h;
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
