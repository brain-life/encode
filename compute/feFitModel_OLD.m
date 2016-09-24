function [fit, w, R2] = feFitModel_OLD(M,dSig,fitMethod)
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
%  - 'lsqnonneg' - MatLab defaoult non-negative least-square solver (SLOW)
%  - 'sgd', 'sgdnn' - Stochastic gradient descent.
%  - 'sgdl1','sgdl1nn' - Stochastic gradient descent with L1 constrain on weights.
%
% See also: feCreate.m, feConnectomeBuildModel.m, feGet.m, feSet.m
%
% Example:
%
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
switch fitMethod
   case {'bbnnls'}
    [nFibers] = size(M,2); %feGet(fe,'nfibers');
    %[nTheta]  = size(M.DictSig,1);
    %[nAtoms] = size(M.DictSig,2); %feGet(fe,'natoms');
    %[Nvoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');
    
    tic
    fprintf('\nLiFE: Computing least-square minimization with BBNNLS...\n')
    opt = solopt;
    opt.maxit = 5000;
    opt.use_tolo = 1;
    out_data = bbnnls_OLD(M,dSig,zeros(nFibers,1),opt);
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
%   case {'sgd','sgdnn'}% stochastic gradient descend, or non-negative stochastic gradient descend
%     tic
%     % Stochastic gradient descent method.
%     % it solves an L2 minimization problem with non-negative constrain.
%     %
%     % Basically it takes 'chuncks' of rows of the M matrix and solves those
%     % separately but contraining to obtain a consistent global solution.
%     signalSiz = size(M,1);
%     if signalSiz >= 1000000
%       siz     = floor(signalSiz * .1); % size of the chuncks (number rows) taken at every iteration of the solver
%     elseif signalSiz > 10000 || signalSiz < 1000000
%       siz     = floor(signalSiz * .5); % size of the chunks (number rows) taken at every iteration of the solver
%     elseif signalSiz <= 10000
%       siz     = signalSiz; % size of the chuncks (number rows) taken at every iteration of the solver
%     else
%       keyboard
%     end
%     stepSiz      = 0.0124; % step in the direction of the gradient, the larger the more prone to local minima
%     stopCriteria = [.1 5 1]; % Stop signals:
%     % First, if total error has not decreased less than
%     %        an XXX proportion of XXXX.
%     % Second, number of small partial fits before
%     %         evaluating the quality of the large fit.
%     % Third, Amount of R2 improvement judged to be
%     %        useful.
%     %        It used to be:  percent improvement in R2
%     %        that is considered a change in quality
%     %        of fit, e.g., 1=1%.
%     n      = 100;       % Number of iteration after which to check for total error.
%     nonneg = strcmpi(fitMethod(end-2:end),'dnn');
%     fprintf('\nLiFE: Computing least-square minimization with Stochastic Gradient Descent...\n')
%     [w, R2] = sgd(dSig,M,siz,        stepSiz,      stopCriteria,        n,         nonneg);
%              %sgd(y,   X,numtoselect,finalstepsize,convergencecriterion,checkerror,nonneg,alpha,lambda)
%     % Save out the Stochastic Gradient Descent parameters
%     fit.params.stepSiz      = stepSiz;
%     fit.params.stopCriteria = stopCriteria;
%     fit.params.numInters    = n;
%     
%     % Save the state of the random generator so that the stochasit cfit can be recomputed.
%     defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
%     fit.randState = defaultStream.State;   
%     
%     % Save out some results 
%     fit.results.R2        = R2;
%     fit.results.nParams   = size(M,2);
%     fit.results.nMeasures = size(M,1);
%     fprintf(' ...fit process completed in %2.3fs\n',toc)
% 
%     case {'sgdl1','sgdl1nn'}% stochastic gradient descend, or non-negative stochastic gradient descend
%     tic
%     % Stochastic gradient descent method.
%     % it solves an L2 minimization problem with non-negative constrain.
%     %
%     % Basically it takes 'chuncks' of rows of the M matrix and solves those
%     % separately but contraining to obtain a consistent global solution.
%     signalSiz = size(M,1);
%     if signalSiz >= 1000000
%       siz     = floor(signalSiz * .1); % size of the chunks (number rows) taken at every iteration of the solver
%     elseif signalSiz > 10000 || signalSiz < 1000000
%       siz     = floor(signalSiz * .5); % size of the chunks (number rows) taken at every iteration of the solver
%     elseif signalSiz <= 10000
%       siz     = signalSiz; % size of the chunks (number rows) taken at every iteration of the solver
%     else
%       keyboard
%     end
%     stepSiz      = 0.0124; % step in the direction of the gradient, the larger the more prone to local minima
%     stopCriteria = [.1 5 1]; % Stop signals:
%     % First, if total error has not decreased less than
%     %        an XXX proportion of XXXX.
%     % Second, number of small partial fits before
%     %         evaluating the quality of the large fit.
%     % Third, Amount of R2 improvement judged to be
%     %        useful.
%     %        It used to be:  percent improvement in R2
%     %        that is considered a change in quality
%     %        of fit, e.g., 1=1%.
%     n      = 100;       % Number of iteration after which to check for total error.
%     nonneg = 1;
%     fprintf('\nLiFE: Computing least-square minimization (L1) with Stochastic Gradient Descent...\n')
%     %lambda = [length(dSig)*2.75];
%     [w, R2] = sgdL1(dSig,M,siz, stepSiz, stopCriteria, n,nonneg,[],lambda);
%     fprintf('Lambda: %2.2f | nFibers: %i | L1 penalty: %2.3f | L2 penalty: %2.3f\n',lambda, length(find(w>0)),sum(w),sum(w.^2))
% 
%     % Save out the Stochastic Gradient Descent parameters
%     fit.params.stepSiz      = stepSiz;
%     fit.params.stopCriteria = stopCriteria;
%     fit.params.numInters    = n;
%     
%     % Save the state of the random generator so that the stochasit cfit can be recomputed.
%     defaultStream = RandStream.getGlobalStream; %RandStream.getDefaultStream;
%     fit.randState = defaultStream.State;   
%     
%     % Save out some results 
%     fit.results.R2        = R2;
%     fit.results.nParams   = size(M,2);
%     fit.results.nMeasures = size(M,1); 
%     fit.results.l2        = sum(w.^2);
%     fit.results.l1        = sum(w);
%     
%     fprintf(' ...fit process completed in %2.3fs\n',toc)
% 
   otherwise
     error('Cannot fit LiFE model using method: %s.\n',fitMethod);
end

% Save output structure.
fit.weights             = w;
fit.params.fitMethod    = fitMethod;

end
