function [fe, results] = FitFullModel(dwiFile, fgFileName, feFileName, L, p, alpha_v, alpha_f, lambda_a, lambda_r)
% INPUT
% dwFile:               diffusion measurements
% fgFileName:           Tractography file
% L:                    discretization parameter in ENCODE
% p:                    Training set ratio, (1-p) is validation set ratio. We divide gradient
%                       directions is training - validation sets to tune parameter lambda
% alpha_v:              l2 regularization parameter single voxel fit
% alpha_f:              l2 regularization parameter full connectome fitting
% lambda_a:             axial diffusivity (lambda1)
% lambda_r:             radial diffusivity (lambda2=lambda3)

%% Initialize the model
tic
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[] ,[], [], L, [lambda_a,lambda_r],0); % We set dwiFileRepeat =  run 02
disp(' ')
disp(['Time for model construction ','(L=',num2str(L),')=',num2str(toc),'secs']);

%% Define training-validation directions by random
nTheta = feGet(fe,'nbvals');
[nAtoms] = feGet(fe,'natoms');
[nFibers] = feGet(fe,'nfibers');
[nVoxels] = feGet(fe,'nvoxels');

%nVoxels = 10000; % For testing purposes ONLY

ind_dirs = randperm(nTheta); 

nTrain = round(p*nTheta); % Number of training directions
nVal = nTheta - nTrain; % Number of validation directions

ind_train = ind_dirs(1:nTrain); % Set of training directions
ind_val = ind_dirs(nTrain+1:end); % Set of validation directions

%% Fit model to TRAINING DATA
Y = fe.life.diffusion_signal_img(1:nVoxels,ind_train)';
D = fe.life.M.Dict(ind_train,:);
Phi = fe.life.M.Phi(:,1:nVoxels,:);


%% STAGE 1: Voxelwise fitting. Fit B and s0 to measurements (alternate between B and s0, see algorithm in Methods_and_Supp.pdf)
% We stop iterating when the maximum number of iterations (Niter) is
% reached or when the error is below a threshold.
Niter = 50;
threshold = 1e-8;

Error_vs_iter = zeros(1,Niter);

B = ttv(Phi,ones(nFibers,1),3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
normY = norm(Y,'fro');
s0 = zeros(nVoxels,1);

Error = norm(Y - ones(nTrain,1)*s0' - D*B,'fro')/normY;

delta = Inf;
n = 1;
disp(' ')
while (n<= Niter)&&(delta > threshold)
    disp(['IN LOOP iter ',num2str(n),' Error=', num2str(Error), ' nnz(B)=',num2str(nnz(B)/numel(B)),' delta=',num2str(delta) ])
    %fprintf('.');
    
    % Min over B
    B = Min_over_B(Y - ones(nTrain,1)*s0',B,D,alpha_v);
    
    % Min over s0
    s0 = Min_over_s0(Y - D*B);
    
    Error_vs_iter(n) = Error;
    Error_ant = Error;
    Error = norm(Y - ones(nTrain,1)*s0' - D*B,'fro')/normY;
    delta = abs(Error - Error_ant);
    n = n + 1;
end

% Reconstruct Phi from B
C = double(sptenmat(ttv(Phi,ones(nAtoms,1),1),[1])); % sum over atoms

b0 = (ones(1,nAtoms)*B)';


S0 = mean(fe.life.diffusion_S0_img(1:nVoxels,:),2);
[i,j,~] = find(C);
values = S0(i);
A0 = sparse(i, j, values, size(C,1), size(C,2));



A = [A0; alpha_f*speye(size(C,2))];

b = [b0; sparse(size(C,2),1)];

opt = solopt;
opt.maxit = 10000; % increase number of iteration to asure convergence
opt.verbose =1;
tic
out = bbnnls_orig(A, b, zeros(nFibers,1), opt);
fprintf('Global Fitting took: %2.3fs.\n',toc)

w = out.x;

error_b = norm(b - A*w)/norm(b);
disp(' ');
disp([' Error_B=', num2str(error_b),' nnz(w)=', num2str(nnz(w))]);

%Phi = reconstruct_Phi(Phi,S0);
B = ttv(Phi,w,3);
[ind, val] = find(B);
B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);

%% Compute Errors and save results
Ypred_train = ones(nTrain,1)*s0' + D*B;
Error_train = norm(Y - Ypred_train,'fro')/normY;
rmse_train = sqrt(mean((Y - Ypred_train).^2,1))./(mean(fe.life.diffusion_S0_img(1:nVoxels,:),2)');
SSres_train = (norm(Y - Ypred_train,'fro'))^2;
Ym = mean(Y(:));
SStot_train = (norm(Y - Ym*ones(size(Y)),'fro'))^2;
R2_train = 1 - SSres_train/SStot_train;

%% Compute error on VALIDATION DATA
Y_val = fe.life.diffusion_signal_img(1:nVoxels,ind_val)';
D_val = fe.life.M.Dict(ind_val,:);
Ypred_val = ones(nVal,1)*s0' + D_val*B;
Error_val = norm(Y_val - Ypred_val,'fro')/norm(Y_val,'fro');
rmse_val = sqrt(mean((Y_val - Ypred_val).^2,1))./(mean(fe.life.diffusion_S0_img(1:nVoxels,:),2))';
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
results.lambda_a = lambda_a;
results.lambda_r = lambda_r;
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


function [Phi] = reconstruct_Phi(Phi, S0)
[nAtoms] = size(Phi,1);
[sub, ~] = find(Phi);
Count = sptensor(sub, ones(size(sub,1),1), size(Phi));
Count = ttv(Count,ones(nAtoms,1),1); % sum entries over mode 1
%[subCount, val] = find(Count);
%Count = sparse(subCount(:,1),subCount(:,2),val,nVoxels,nFibers); % Sparse matrix (nVoxels x nFibers) containing the number of atoms per voxel per fiber

Count = double(sptenmat(Count,[1])); % Sparse matrix (nVoxels x nFibers) containing the number of atoms per voxel per fiber



ind = sub2ind(size(Count), sub(:,2), sub(:,3));
count = Count(:);
%count = double(sptenmat(Count,[1,2])); % equivalent to count = Count(:), sparse vectorization
div = full(count(ind));

%valPhi = w(sub(:,3)).*S0(sub(:,2)); % assign weight to fascicle slice
valPhi = S0(sub(:,2)); % assign weight to fascicle slice
valPhi = valPhi./div;


Phi = sptensor(sub, valPhi, size(Phi));

end

 
function [B] = Min_over_B(Y,B,D,lambda)
nVoxels = size(Y,2); % Number of voxels
Bnew = zeros(size(B));
parfor v=1:nVoxels % we use parfor because each voxel is independent of the rest of voxels
    [ind, val] = find(B(:,v)); % each instance of par for works on a column in matrix B
    
    C = [D(:,ind); lambda*eye(length(ind))]; % augmented matrix for Tikhonov regularizer
    d = [Y(:,v); zeros(length(ind),1)]; % augmented vector
    
    % Set parameters for the BBNLS algorithm
    opt = solopt;
    out = bbnnls_orig(C, d, zeros(size(C,2),1), opt);
    b = out.x;
    
    col = zeros(size(B,1),1);
    col(ind,1) = b
    
    Bnew(:,v) = col;
    
end
% Finally, we return a sparse matrix
B = sparse(Bnew);
end


function [s0] = Min_over_s0(E)
s0 = (sum(E,1)/size(E,1))';

s0(s0<=0) = 0;

end

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
