function fe = feConnectomeEncoding(fe)
% Compute multiway decompositon model to predict directional diffusion in each voxel from fibers
%
%   fe = feConnectomeBuildModel(fe)
%
% INPUTS: fe -  An fe structure, see feCreate.m
%
% See also: feFitModel.m, feComputePredictedSignal.m
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%

if notDefined('fe'),  error('LiFE (fe = feCreate) struct needed'); end
if ~isfield(fe,'life')
  error('LiFE - the field ''life'' is necessary in the fe structure.')
end

fprintf('\n[%s] Encoding connectome (building Phi, sparse tensor) ... ',mfilename); 
tic

% make sure that Matlab accept Very Large matrices
s = Settings();
set(s.matlab.desktop.workspace, 'ArraySizeLimitEnabled',false)

nFibers      = feGet(fe,'n fibers');
nTotalNodes = fefgGet(fe.fg,'n total nodes');
fibers = fe.fg.fibers;
imgsize = feGet(fe,'volumesize');
imgsize = imgsize(1:3); % 4th dimension is discarded
nTotalVoxels = prod(imgsize); % including voxels not in the ROI

nAtoms       = feGet(fe,'n Atoms');
%orient       = feGet(fe,'orient');
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

fe = feSet(fe,'Indication Tensor',Phi);

clear 'A' 'vox' 'fib' 'vals'

fprintf('took: %2.3fs.\n',toc)

return
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
