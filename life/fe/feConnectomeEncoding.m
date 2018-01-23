function fe = feConnectomeEncoding(fe)
% Compute multiway decompositon model to predict directional diffusion in each voxel from fibers
%
%   fe = feConnectomeBuildModel(fe)
%
% INPUTS: fe -  An fe structure, see feCreate.m
%
% See also: feFitModel.m, feComputePredictedSignal.m
%
%  Copyright (2017), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%

if notDefined('fe'),  error('LiFE (fe = feCreate) struct needed'); end
if ~isfield(fe,'life')
  error('LiFE - the field ''life'' is necessary in the fe structure.')
end

fprintf('\n[%s] Encoding connectome (building Phi, sparse tensor) ... ',mfilename); 
tic

% make sure that Matlab accept Very Large matrices
if verLessThan('matlab','9.2')
    % -- Code to run in MATLAB R2016b and earlier here --
    s = Settings();
    set(s.matlab.desktop.workspace, 'ArraySizeLimitEnabled',false);
end

% Set MAXMEM available
MAXMEM = getenv('MAXMEM'); % read setting from enviroment in kb
if MAXMEM
    disp(['MAXMEM set to ',MAXMEM]);
    MAXMEM = str2num(MAXMEM);
elseif (isunix||ismac)
    disp('MAXMEM not set, need to calculate (UNIX or MacOS)')
    [~,out] = system('cat /proc/meminfo |grep MemFree');
    textcell = regexp(out,'\d*','Match');
    MAXMEM = str2num(textcell{1});
elseif ispc
    disp('MAXMEM not set, need to calculate (PC)')
    user = memory;
    MAXMEM = user.MaxPossibleArrayBytes/1024; % Max mem for arrays in Kb
end

if ~MAXMEM 
    error('Could not determine MAXMEM');
end


% Check required number of nodes and split the tensor Phi computation in
% pieces having max nNodesMax nodes each

%nNodesMax = 30000000; % Maximum nodes per batch. This is a reference value, for example one HCP3T subject has 28,677,744 nodes.
%nNodesMax = 12000000; % For testing

%nNodesMax = round(MAXMEM*30000000/32000000); % In Karst we have available 32,000,000kb of memory, which allows to process up to 30,000,000 nodes;

nNodesMax = round(MAXMEM*15000000/32000000); % 

nTotalNodes = fefgGet(fe.fg,'n total nodes');
nFibers      = feGet(fe,'n fibers');
nBatch = ceil(nTotalNodes/nNodesMax); % Number of batch
nFib_Batch = ceil(nFibers/nBatch); % Number of fibers per batch

disp(['nNodesMax =',num2str(nNodesMax),', Total number of nodes = ',num2str(nTotalNodes),',  number of batch computation = ',num2str(nBatch)])

fprintf('\n');
for n=1:nBatch
    fprintf('Encoding batch % 2.f\n',n)
    fibers_range = (n-1)*nFib_Batch + 1: min(n*nFib_Batch,nFibers);
    if n == 1
        Phi = compute_Phi_batch(fe,fibers_range); % Phi is created in the first batch
    else
        Phi = concatenate_mode3(Phi, compute_Phi_batch(fe,fibers_range));
    end
end

fe = feSet(fe,'Indication Tensor',Phi);

fprintf('took: %2.3fs.\n',toc)

return
end

function [Phi1] = concatenate_mode3(Phi1, Phi2)
subs1 = Phi1.subs;
vals1 = Phi1.vals;
s1 = size(Phi1);
s2 = size(Phi2);

Nvals1 = length(vals1);
subs1 = [subs1; Phi2.subs];
vals1 = [vals1; Phi2.vals];

subs1(Nvals1+1:end,3) = subs1(Nvals1+1:end,3) + s1(3);
s1(3) = s1(3) + s2(3);
Phi1 = sptensor(subs1,vals1,s1);

end

function [Phi] = compute_Phi_batch(fe,fibers_range)
nFibers      = length(fibers_range);
fg = fe.fg;
fg.fibers = fg.fibers(fibers_range);
nTotalNodes = fefgGet(fg,'n total nodes');
fibers = fe.fg.fibers(fibers_range);
imgsize = feGet(fe,'volumesize');
imgsize = imgsize(1:3); % 4th dimension is discarded
nTotalVoxels = prod(imgsize); % including voxels not in the ROI

nAtoms       = feGet(fe,'n Atoms');
Nphi         = feGet(fe,'Nphi');
Ntheta       = feGet(fe,'Ntheta');

% Compute fiber
[tubes,grad] = fiberOfNodes(fibers,nTotalNodes); % this function gives a (1xnTotalNodes) containing the fiber number for each node
fibers = cell2mat(fibers(:)'); 

% Compute voxels
voxel_coord = round(fibers) + 1;
cols = sub2ind(imgsize, voxel_coord(1,:)', voxel_coord(2,:)', voxel_coord(3,:)');

% Compute atoms
%grad = gradient(fibers);
grad = grad./repmat(sqrt(sum(grad.^2)),3,1); % normalize columns

rows = get_atom(grad,Nphi,Ntheta); 

% Construct Indication Tensor Phi
Phi = sptensor([rows,cols,tubes],ones(nTotalNodes,1),[nAtoms,nTotalVoxels,nFibers]);

% The following sparse matrix (Nvoxels x Nfibers) counts the number of
% nodes per fiber per voxel.
A = sparse(cols,tubes,ones(nTotalNodes,1),nTotalVoxels,nFibers);

% Get voxex indices in roi
roi_coords = feGet(fe,'roicoords');
roi_ind = sub2ind(imgsize,roi_coords(:,1)',roi_coords(:,2)',roi_coords(:,3)');

% Restrict the Sparse Tensor and matrix A to the ROI voxels only
Phi = Phi(:,roi_ind,:); % reduce tensor in 3rd dimension to roi voxels only
A = A(roi_ind,:);

vox = Phi.subs(:,2); % Recover voxel indices
fib = Phi.subs(:,3); % Recover fiber indices
vals = Phi.vals; % Recover vals

% Normalize values in the sparse tensor dividing by the number of nodes per
% voxel and per fiber
nVoxels = length(roi_ind);
a = A(:);
vals = vals./a(sub2ind([nVoxels,nFibers],vox,fib));

% Create the Sparse Tensor (normalized)
Phi = sptensor(Phi.subs,vals,size(Phi));

% The following multiplies every slice by the corresponding S0(voxel) value
S0 = feGet(fe,'s0_img');
vals = Phi.vals.*S0(Phi.subs(:,2));
Phi = sptensor(Phi.subs,vals,size(Phi));

end




function [tubes, grad] = fiberOfNodes(fibers,nTotalNodes)
    tubes = zeros(nTotalNodes,1);
    grad = zeros(3,nTotalNodes);
    node = 1;
    for f = 1:size(fibers,1)
        NumberOfNodes = size(fibers{f},2);
        tubes(node:node+NumberOfNodes-1) = f;
        grad(:,node:node+NumberOfNodes-1) = gradient(fibers{f});
        node = node + NumberOfNodes;
    end
end


function atom_ind = get_atom(vect,Nphi,Ntheta)
    deltaTheta = pi/Ntheta;
    deltaPhi = pi/Nphi;
    
    % vector should belong to the positive half of the sphere
    vect(:,vect(2,:)<0) = -vect(:,vect(2,:)<0);
    
    [angPhi,angTheta,r] = cart2sph(vect(1,:)',vect(2,:)',vect(3,:)');
    angTheta = pi/2 - angTheta; % we measure theta as the angle with the positive semi-axis y
    % ind_Phi=1 correspond to Phi=0, ind_Phi=Nphi to Phi = Pi - delta_Phi
    ind_Phi = round(angPhi/deltaPhi) + 1;
    
    % ind_Phi=1 correspond to Phi=0, ind_Phi=Nphi to Phi = Pi - delta_Phi
    indsub = find(ind_Phi == Nphi+1);
    ind_Phi(indsub) = 1;
    angTheta(indsub) = pi - angTheta(indsub);
    
    % ind_Theta=1 correspond to Theta=delta_Theta, ind_Theta=NTheta-1 to
    % Phi = Pi - delta_Theta
    ind_Theta = round(angTheta/deltaTheta);
    
    indsub = find((ind_Theta ~= 0) & (ind_Theta ~= Ntheta));
    atom_ind = ones(size(ind_Phi));
    atom_ind(indsub) = sub2ind([Ntheta-1,Nphi],ind_Theta(indsub),ind_Phi(indsub)) + 1;
    
    
end
