
function [cross_vox_ind, cross_vox_coord] = feFindCrossingVoxCoord(fe, ind1, ind2)

[Na] = size(fe.life.M.Phi,1); % # of atoms
[Nv] = size(fe.life.M.Phi,2); % # of voxels
[Nf] = size(fe.life.M.Phi,3); % # of fascicles

Phi_tract1 = fe.life.M.Phi(:,:,ind1);
[subs, vals] = find(Phi_tract1);
vox_ind_1 = unique(subs(:,2));

Phi_tract2 = fe.life.M.Phi(:,:,ind2);
[subs, vals] = find(Phi_tract2);
vox_ind_2 = unique(subs(:,2));

cross_vox_ind = intersect(vox_ind_1, vox_ind_2);
cross_vox_coord = fe.roi.coords(cross_vox_ind,:);


end
