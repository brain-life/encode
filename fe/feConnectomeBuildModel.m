function fe = feConnectomeBuildModel(fe,Compute_matrix_M)
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

fprintf('\n[%s] Building the Connetome Model (Indication Tensor Phi) ... ',mfilename); 
tic


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
voxel_coord = ceil(fibers) + 1;
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

roi_coords = feGet(fe,'roicoords');
roi_ind = sub2ind(imgsize,roi_coords(:,1)',roi_coords(:,2)',roi_coords(:,3)');

Phi = Phi(:,roi_ind,:); % reduce tensor in 3rd dimension to roi voxels only
A = A(roi_ind,:);

vox = Phi.subs(:,2);
fib = Phi.subs(:,3);
vals = Phi.vals;

nVoxels = length(roi_ind);
a = A(:);
vals = vals./a(sub2ind([nVoxels,nFibers],vox,fib));

Phi = sptensor(Phi.subs,vals,size(Phi));

% The following multiplies every slice by the corresponding S0(voxel) value
S0 = feGet(fe,'s0_img');
vals = Phi.vals.*S0(Phi.subs(:,2));
Phi = sptensor(Phi.subs,vals,size(Phi));

fe = feSet(fe,'Indication Tensor',Phi);


if Compute_matrix_M % This was introduced in order to compare OLD vs NEW LiFE
    % Compute Large sparse matrix M (as in the old LiFE)
    disp('Computing large and sparse matrix M (as in the old LiFE)')
    
    nBvecs       = feGet(fe,'nBvecs');
    bvecs        = feGet(fe,'bvecs');  % bvecs
    bvals        = feGet(fe,'bvals');  % bvals
    nodeSig = zeros(nBvecs,nTotalNodes);

    D = diag(fe.life.modelTensor); % diagonal matix with diffusivities
    parfor j=1:nTotalNodes
        disp(['Computing diffusion of node ',num2str(j),'/',num2str(nTotalNodes)]);
        [Rot,~, ~] = svd(grad(:,j)); % Compute the eigen vectors of the kernel orientation 
        Q = Rot*D*Rot';
        nodeSig(:,j) = exp(- bvals .* diag(bvecs*Q*bvecs')); % Compute the signal contribution of a fiber in the kernel orientation divided S0
        nodeSig(:,j) = nodeSig(:,j) - mean(nodeSig(:,j)); % demeaned signal
    end

    indi = zeros(nBvecs*nTotalNodes,1);
    indj = zeros(nBvecs*nTotalNodes,1);
    vals = zeros(nBvecs*nTotalNodes,1);

    % Construction of nonzero indices of matrix M    
    for node=1:nTotalNodes
        disp(['Building matrix M, node ',num2str(node),'/',num2str(nTotalNodes)]);
        indi((node-1)*nBvecs + 1 : node*nBvecs) = [(cols(node)-1)*nBvecs+1:cols(node)*nBvecs]';
        indj((node-1)*nBvecs + 1 : node*nBvecs) = repmat(tubes(node),[nBvecs,1]);
        vals((node-1)*nBvecs + 1 : node*nBvecs) = nodeSig(:,node);   
    end
    
    % June 27/6/2015
    clear nodeSig Phi rows cols tubes fibers grad voxel_coord

    Mmatrix = sparse(indi, indj, vals);
    kept_ind = zeros(length(roi_ind)*nBvecs,1);
    for i=1:length(roi_ind)
        disp(['Keeping voxels that are in ROI, voxel= ',num2str(i),'/',num2str(length(roi_ind))]);
        kept_ind((i-1)*nBvecs+1:i*nBvecs) = (roi_ind(i)-1)*nBvecs + 1: roi_ind(i)*nBvecs;
    end

    Mmatrix = Mmatrix(kept_ind,:); % reduce Matrix
    [indi,indj,vals]=find(Mmatrix);

    vox = ceil(indi/nBvecs);
    vals = vals.*S0(vox);
    Mmatrix = sparse(indi,indj,vals);


    fe.life.M.Mmatrix = Mmatrix;
end

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


% This function computes the atom indexing in a (Nphi x Ntheta)-grid dictionary corresponding to a given
% orientation vect
% function atom_ind = get_atom(vect,Nphi,Ntheta)
%     deltaTheta = pi/Ntheta;
%     deltaPhi = pi/Nphi;
%     
%     % vector should belong to the positive half of the sphere
%     if vect(2) < 0
%         vect = - vect;
%     end
%     
%     % cartesian to spherical coordinates
%     [angPhi,angTheta,r] = cart2sph(vect(1),vect(2),vect(3));
%     angTheta = pi/2 - angTheta; % we measure theta as the angle with the positive semi-axis y
%     
%     % ind_Phi=1 correspond to Phi=0, ind_Phi=Nphi to Phi = Pi - delta_Phi
%     ind_Phi = round(angPhi/deltaPhi) + 1;
%     if ind_Phi == Nphi+1
%         ind_Phi = 1;
%         angTheta = pi - angTheta;
%     end
% 
%     % ind_Theta=1 correspond to Theta=delta_Theta, ind_Theta=NTheta-1 to
%     % Phi = Pi - delta_Theta
%     ind_Theta = round(angTheta/deltaTheta);
%      
%     if (ind_Theta == 0) || (ind_Theta == Ntheta)
%         atom_ind = 1; % spin-Up or spin-Down case
%     else
%         atom_ind = sub2ind([Ntheta-1,Nphi],ind_Theta,ind_Phi) + 1;
%     end
% end
