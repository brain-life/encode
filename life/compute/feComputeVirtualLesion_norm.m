% This function compute the rmse in a path neighborhood voxels with and
% without Virtual Lesion
function [ rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = feComputeVirtualLesion_norm(fe, ind_tract)
% INPUTS:
% fe: fe structure
% ind1: indices to fibers in the tract to be virtually lesioned

% ind_nnz = find(fe.life.fit.weights);
% ind_tract = ind_nnz(ind1);

% We want to find which voxels that a group of fibers crosses.
[inds, vals] = find(fe.life.M.Phi(:,:,ind_tract)); % find nnz entries of subtensor

% inds has a list of (i,j,k) positions of nnz entries. Since we are interested in
% locating voxels we need to look at the second column (j).

voxel_ind = unique(inds(:,2));

% To find which other fibers crosses these voxels, we need to look at
% the the subtensor that corresponds to those voxels
% See following lines


ind2 = feGet(fe,'pathneighborhood',ind_tract);
nFib_tract = length(ind_tract);
 
nVoxels = length(voxel_ind);
nTheta = feGet(fe,'n bvals');

% [inds, vals] = find(fe.life.M.Phi(:,voxel_ind,:)); % find indices for nnz in the subtensor defined by voxel_ind
% ind2 = unique(inds(:,3)); % Fibers are the 3rd dimension in the subtensor

w = fe.life.fit.weights;
ind_nnz = find(w);
nFib_PN = length(ind2);

measured  = feGet(fe,'dsigdemeaned by voxel');

%% DEBUGGING Plot path, path-neightboord and voxels
% figure('name','path-nehighborhood and tract','color','w')
% hold on
% fg = fe.fg.fibers;
% for ii = 1:length(ind_tract)
%     plot3(fg{ind_tract(ii)}(1,:),fg{ind_tract(ii)}(2,:),fg{ind_tract(ii)}(3,:),'r-');
% end
% view(90,0)

% ind_nnz = find(w);
% %tmp_ind = randsample(ind2,200);
% tmp_ind = ind2;
% tmp_ind = intersect(tmp_ind,ind_nnz);
% for ii = 1:length(tmp_ind)
%  plot3(fg{tmp_ind(ii)}(1,:),fg{tmp_ind(ii)}(2,:),fg{tmp_ind(ii)}(3,:),'b-')    
% end
% drawnow;

% c = fe.roi.coords - 1.5*ones(size(fe.roi.coords));
% c = c(voxel_ind,:);
% %tmp_c_ind = randsample(1:size(c,1),200);
% tmp_c_ind = 1:size(c,1);
% plot3(c(tmp_c_ind,1),c(tmp_c_ind,2),c(tmp_c_ind,3),'co')

%% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
S0 = fe.life.diffusion_S0_img(voxel_ind);
measured = measured(:,voxel_ind);
% Restrict tensor model to the PN voxels
M = fe.life.M;
M.Phi = M.Phi(:,voxel_ind,:);
predicted_woVL =  reshape(M_times_w(M.Phi.subs(:,1),M.Phi.subs(:,2),M.Phi.subs(:,3),M.Phi.vals,M.DictSig,w,nTheta,nVoxels),size(measured));
%predicted_woVL =  reshape(M_times_w(M,w),size(measured));
rmse_woVL = sqrt(mean((measured - predicted_woVL).^2,1));
rmse_woVL = rmse_woVL./S0';

%% Compute rmse restricted to Path neighborhood voxels with Virtual Lesion
w_VL = w;
w_VL(ind_tract) = 0;
predicted_VL =  reshape(M_times_w(M.Phi.subs(:,1),M.Phi.subs(:,2),M.Phi.subs(:,3),M.Phi.vals,M.DictSig,w_VL,nTheta,nVoxels),size(measured));
rmse_wVL = sqrt(mean((measured - predicted_VL).^2,1));
rmse_wVL = rmse_wVL./S0';


end
