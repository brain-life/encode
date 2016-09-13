function fe = feGetConnectomeInfo(fe)
% Find the unique fibers in a connectome.
% 
%   unique_fibers_index = feGetConnectomeInfo(fe)
%
% We identify the fibers that have identical trajectories (that pass
% through exactly the same voxels) within the ROI. These are stored in the
% fe slot 'unique index'.
%
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Get the indexes to the voxels actually used to build LifE.
% Remember that voxels with no fibers are disregarded during the build.
usedVoxels = feGet(fe,'usedVoxels');
nVoxels    = length(usedVoxels);

% preallocate memory for speed.
f                  = cell( nVoxels,1); % This will contain all the fibers in each voxel
unique_fibers_index= cell( nVoxels,1); % This will contain the *unique* fibers in each voxel
unique_fibers_num  = zeros(nVoxels,1); % This will contain the number of unique fibers in each voxel

% Get the fibers and unique fibers from voxel2FNpairs, in each voxel:
parfor vv = 1:nVoxels
  % This is the index of the voxel used, it is used to address voxel2FNpair
  voxIndex = usedVoxels(vv);
  
  % Indexes to the fibers in the current voxel
  f{vv} = fe.life.voxel2FNpair{voxIndex}(:,1);
  
  % Reduce to the unique fibers:
  unique_fibers_index{vv} = sort(unique(f{vv}));
  unique_fibers_num(vv)   = length(unique_fibers_index{vv});
end

tot_fibers_num = cellfun(@length,f);

% Set them in the fe structure: 
% Indexes to the unique fibers going through each voxel
fe = feSet(fe,'index to unique fibers in each voxel',  unique_fibers_index);
% Number of unique fibers going through each voxel
fe = feSet(fe,'number of unique fibers in each voxel', unique_fibers_num);

% Number of total fibers going trhough each voxel
fe = feSet(fe,'number of total fibers in each voxel',  tot_fibers_num);
% Indexes to the total fibers going through each voxel
fe = feSet(fe,'index of total fibers in each voxel',   f);

return
