function [fe, results] = FitFullModelSampleAllTracts(dwiFile, fgFileName, feFileName, L, p, n, alpha_v, alpha_f, lambda_1, lambda_2, fg_classification)
% INPUT
% dwFile:               diffusion measurements
% fgFileName:           Tractography file
% L:                    discretization parameter in ENCODE
% p:                    Training set ratio, (1-p) is validation set ratio. We divide gradient
%                       directions is training - validation sets to tune parameter lambda
% n:                    % of voxels randomly subsample for fitting
% alpha_v:              l2 regularization parameter single voxel fit
% alpha_f:              l2 regularization parameter full connectome fitting
% lambda_1:             axial diffusivity (lambda1)
% lambda_2:             radial diffusivity (lambda2=lambda3)

%% Initialize the model
tic
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[] ,[], [], L, [lambda_1,lambda_2],0); % We set dwiFileRepeat =  run 02
disp(' ')
disp(['Time for model construction ','(L=',num2str(L),')=',num2str(toc),'secs']);

%% Define training-validation directions by random
nTheta = feGet(fe,'nbvals'); % Number of gradient directions
ind_dirs = randperm(nTheta); % Random permutation of gradient directions indices

nTrain = round(p*nTheta); % Number of training directions
nVal = nTheta - nTrain; % Number of validation directions

ind_train = ind_dirs(1:nTrain); % Set of training directions
ind_val = ind_dirs(nTrain+1:end); % Set of validation directions

nAtoms = feGet(fe,'natoms'); % Number of atoms in the dictionary
nFibers = feGet(fe,'nfibers'); % Number of Fibers in the Connectome

%% Find voxels for All Tracts
load(fg_classification); % load classification of fascicles
ind_fibers_in_Tracts = find(classification.index~=0); % find all fibers classified as belonging to a tract
ind_vox_tracts = feGet(fe,'vox ind from fibers',ind_fibers_in_Tracts); % get voxel indices for those fascicles
nVoxels = length(ind_vox_tracts);

%% Subsampling of voxels
nVoxSample = round(n*nVoxels); % we subsample within the voxels of tracts
ind_vox_within_tracts = randperm(nVoxels,nVoxSample); % random selection of few voxels

% Determine indices of subsampled voxels
ind_vox = ind_vox_tracts(ind_vox_within_tracts);
disp(['Sampled voxels: ', num2str(100*nVoxSample/feGet(fe,'nvoxels')),'%, ', num2str(nVoxSample), ' voxels'])

%% Fit model to TRAINING DATA
% First, we define our variables of ENCODE model (see eq. (3), Methods_and_Supp.pdf) restricted to selected directions and voxels:
Y = fe.life.diffusion_signal_img(ind_vox,ind_train)'; % diffusion data matrix restricted to selected voxels and directions
normY = norm(Y,'fro'); % Frobenius norm of diffusion data matrix Y
D = fe.life.M.Dict(ind_train,:); % Dictionary of diffusion kernels restricted to selected directions
Phi = fe.life.M.Phi(:,ind_vox,:); % Sparse Core tensor restricted to the selected voxels

%% STAGE 1: Voxelwise fitting. Fit B and s0 to measurements (alternate between B and s0, see algorithm in Methods_and_Supp.pdf)
% We stop iterating when the maximum number of iterations (Niter) is
% reached or when the error is below a threshold.
Niter = 50; % Number of iterations
threshold = 1e-8; % Error threshold

% Initialize B using eq. (6) in Methods_and_Supp.pdf, assuming an all-ones
% weights vector w.
B = ttv(Phi,ones(nFibers,1),3); 
% We need to convert B from tensor to sparse matrix format
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxSample);

% Initialize s0 with all-zeros (s0 is the isotropic diffusion which is indicated as I_0 in Methods_and_Supp.pdf) 
s0 = zeros(nVoxSample,1);

Error_vs_iter = zeros(1,Niter);

% Compute initial relative Error of the model
Error = norm(Y - ones(nTrain,1)*s0' - D*B,'fro')/normY;

% Main Loop
delta = Inf;
n = 1;
disp(' ');
while (n<= Niter)&&(delta > threshold)
    disp(['IN LOOP iter ',num2str(n),' Error=', num2str(Error), ' nnz(B)=',num2str(nnz(B)/numel(B)),' delta=',num2str(delta) ])
    
    %% Min over B (consider s0 fixed)
    %% THIS IS THE FUNCTION WE NEED TO PARALLELIZE
    B = Min_over_B(Y - ones(nTrain,1)*s0', B, D, alpha_v); % eq. (10) in Methods_and_Supp.pdf
    
    %% Min over s0 (consider B fixed)
    s0 = Min_over_s0(Y - D*B); % eqs. (8) and (9) in Methods_and_Supp.pdf
    
    Error_vs_iter(n) = Error;
    Error_ant = Error;
    Error = norm(Y - ones(nTrain,1)*s0' - D*B,'fro')/normY;
    delta = abs(Error - Error_ant);
    n = n + 1;
end

%% STAGE 2: Global fitting, computing weights w 
% We need to construct a Nonnegative Least squares problem with l2 regularization such that we
% minimize ||b - Ax|| + alpha_f*||alpha_f||^2 (see eq. (14) in
% Methods_and_Supp.pdf), which is equivalent to solve a classical NNLS
% problem with extended matrix A and vector b.

% Need to compute auxiliary 
C = double(sptenmat(ttv(Phi,ones(nAtoms,1),1),[1])); % sum over atoms
b0 = (ones(1,nAtoms)*B)';

S0 = mean(fe.life.diffusion_S0_img(ind_vox,:),2); % Compute S0 image by averaging over multiple values
[i,j,~] = find(C);
values = S0(i);
A0 = sparse(i, j, values, size(C,1), size(C,2));

% Now, we are ready to define the extended matrix and vector and solve the
% NNLS problem
A = [A0; alpha_f*speye(size(C,2))]; % Augmented matrix
b = [b0; sparse(size(C,2),1)]; % Augmented vector (this is referenced as z in Methods_and_Supp.pdf)

% We set parameters for the NNLS optimization algorithm
opt = solopt; % default values
opt.maxit = 10000; % increase number of iteration to asure convergence
opt.verbose =0; % hide results
tic 
% We use the BBNNLS algorithm from Kim et al, "A non-negative least squares". Optimization Methods and Software 28, 1012-1039 (2013)
out = bbnnls_orig(A, b, zeros(nFibers,1), opt);
fprintf('Global Fitting took: %2.3fs.\n',toc)

w = out.x; % we set the value of weights

error_b = norm(b - A*w)/norm(b); % compute the error fitting vector b
disp(' ');
disp([' Error_B=', num2str(error_b),' nnz(w)=', num2str(nnz(w))]);

% Update of matrix B with the fitted weights
B = ttv(Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxSample); % convert it to a sparse matrix

%% Compute Errors and save results
Ypred_train = ones(nTrain,1)*s0' + D*B;
Error_train = norm(Y - Ypred_train,'fro')/normY;
rmse_train = sqrt(mean((Y - Ypred_train).^2,1))./(mean(fe.life.diffusion_S0_img(ind_vox,:),2)');
SSres_train = (norm(Y - Ypred_train,'fro'))^2;
Ym = mean(Y(:));
SStot_train = (norm(Y - Ym*ones(size(Y)),'fro'))^2;
R2_train = 1 - SSres_train/SStot_train;

%% Compute error on VALIDATION DATA
Y_val = fe.life.diffusion_signal_img(ind_vox,ind_val)';
D_val = fe.life.M.Dict(ind_val,:);
Ypred_val = ones(nVal,1)*s0' + D_val*B;
Error_val = norm(Y_val - Ypred_val,'fro')/norm(Y_val,'fro');
rmse_val = sqrt(mean((Y_val - Ypred_val).^2,1))./(mean(fe.life.diffusion_S0_img(ind_vox,:),2))';
SSres_val = (norm(Y_val - Ypred_val,'fro'))^2;
Ym = mean(Y_val(:));
SStot_val = (norm(Y_val - Ym*ones(size(Y_val)),'fro'))^2;
R2_val = 1 - SSres_val/SStot_val;

disp(' ');
disp([' Error_TRAIN=', num2str(Error_train), ' Error_VAL=', num2str(Error_val), ' R2_VAL=', num2str(R2_val)]);
disp([' rmse_TRAIN=', num2str(mean(rmse_train)), ' rmse_VAL=', num2str(mean(rmse_val))]);

%% Save results
results.Error_train = Error_train;
results.Error_val = Error_val;
results.rmse_train = rmse_train;
results.rmse_val = rmse_val;
results.R2_train = R2_train;
results.ind_train = ind_train;
results.nTrain = nTrain;
results.ind_val = ind_val;
results.nVal = nVal;
results.alpha_v = alpha_v;
results.alpha_f = alpha_f;
results.lambda_1 = lambda_1;
results.lambda_2 = lambda_2;
results.R2_val = R2_val;
results.Error_B_vs_iter = Error_vs_iter;
results.nnz = length(unique(fe.life.M.Phi.subs(:,3)));
results.p = p;
results.L = L;

fe.life.M.Phi = Phi;
fe.life.s0 = s0;
fe.life.B = B;
fe.life.fit.results = results;
fe.life.fit.weights = w;

end


%% The function below need to be implemented efficiently in parallel to avoid parfor and use,
%% for example GPU to process each voxel in parallel.
% This function provides an estimation of matrix B given the rest of variables in the system.
% Matrix B is solved column by column, where each column correspond to a particular voxel. 
% For each voxel we need to solve a NNLS optimization problem using an extended matrix and vector in order to incorporate the l2 regurarizer controlled by the parameter alpha_v
function [B] = Min_over_B(Y,B,D,lambda)
nVoxels = size(Y,2); % Number of voxels
for v=1:nVoxels % we use parfor because each voxel is independent of the rest of voxels
    [ind, val] = find(B(:,v)); % each instance of par for works on a column in matrix B
    
    C = [D(:,ind); lambda*eye(length(ind))]; % augmented matrix for Tikhonov regularizer
    d = [Y(:,v); zeros(length(ind),1)]; % augmented vector
    
    % Set parameters for the BBNLS algorithm
    opt = solopt;
    tic; out = bbnnls_orig_gpu(C, d, zeros(size(C,2),1), opt); toc
    b = out.x;
    
    Bvals{v}.ind = ind;
    Bvals{v}.val = b;
    
end
% Finally, we reconstract the sparse matrix B from theirs columns
for v=1:nVoxels
    B(Bvals{v}.ind,v) = Bvals{v}.val;
end

end

%% Explicit computation of the isotropic diffusion from matrix E
function [s0] = Min_over_s0(E)
s0 = (sum(E,1)/size(E,1))';

s0(s0<=0) = 0;

end


%% Options for the optimization algorithm BBNNLS
function options = solopt(varargin)
% SOLOPT  --  Creates a default options structure for BBNNLS
%
% OPTIONS = SOLOPT
%

options.asgui = 0;
options.beta = 0.0498;
options.compute_obj = 1;
% diminishing scalar; beta^0 =  opt.dimbeg
% beta^k = opt.dimbeg / k^opt.dimexp
options.dimexp = .5;
options.dimbeg = 5;
options.maxit = 1000;
options.maxtime = 10;
options.maxnull = 10;
options.max_func_evals = 30;
options.pbb_gradient_norm = 1e-9;
options.sigma = 0.298;
options.step  = 1e-4;
options.tau = 1e-7;             
options.time_limit = 0;
options.tolg = 1e-3;
options.tolx = 1e-8;
options.tolo = 1e-5;
options.truex=0;
options.xt=[];
options.use_kkt = 0;
options.use_tolg = 1;
options.use_tolo = 0;
options.use_tolx = 0;
options.useTwo = 0;
options.verbose = 0;                    % initially
if nargin == 1
  options.variant = varargin{1};
else   % Default
  options.variant = 'SBB';
end
end
