% This function return the indices to the atoms having a particular spatial orientation determined by a main_orient +- offest
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
function [ ind] = feGetVoxels(fe, v0, dv)
% INPUTS:
% fe: fe structure
% v0: [x,y,z] coordinates of central voxel
% dv: vecinity size

% OUTPUT:
% ind: indices to voxels located in a vecinity of v0


Nv = size(fe.roi.coords,1);

d2 = sum(fe.roi.coords.^2,2) - 2*fe.roi.coords*v0' + repmat(norm(v0)^2,Nv,1);
ind = find(d2 < dv^2);

end

