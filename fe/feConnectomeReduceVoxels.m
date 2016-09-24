function [fe, indicesFibersKept] = feConnectomeReduceVoxels(fe,voxelsToKeep, fibersToKeep)
% Select the voxels to keep (extract) in a connectome matrix
%
%   [fe, indicesFibersKept] = feConnectomeReduceVoxels(fe,voxelToKeep)
%
% voxelsToKeep:  Is a binary list of voxels we preserve.
%
% We expand the voxelsToKeep into a binary list of 0's and 1's that in
% which each is expanded by nBvecs.  The 1s are the rows of the M matrix we
% will keep.
%
%  Example:
%
% Copyright (2015-2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% Get the indices to each voxels' signal --DON"T NEED
% vxRows = feGet(fe,'voxelrows',voxelsToKeep);

% Return only the mode and the signal for the voxels we want to keep
if (numel(size(fe.life.M.Phi)) == 2)
   fe.life.M.Phi = fe.life.M.Phi(:,voxelsToKeep);
else
   fe.life.M.Phi = fe.life.M.Phi(:,voxelsToKeep, :);
end

fe.life.diffusion_signal_img  = fe.life.diffusion_signal_img(voxelsToKeep,:);

% Set the new number of voxels, by indexing inside the roi and
% returning it as an ROI.
fe.roi.coords = fe.roi.coords(voxelsToKeep,:);

% Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
fe.life.diffusion_S0_img = fe.life.diffusion_S0_img(voxelsToKeep);

% Now remove signals for the second data set if it was loaded
if isfield(fe,'rep')   
    if ~isempty(fe.rep.diffusion_signal_img)
    % Set the new diffusion signal, the one for only these subset of voxels.
    fe.rep.diffusion_signal_img = fe.rep.diffusion_signal_img(voxelsToKeep,:);
    end
    
    if ~isempty(fe.rep.diffusion_S0_img)
    % Set the diffusion signal at 0 diffusion weighting (B0) for this voxel:
    fe.rep.diffusion_S0_img = fe.rep.diffusion_S0_img(voxelsToKeep); 
    end
    
end

% Remove the fibers of the fascicle from the fe.
fe = feConnectomeReduceFibers(fe, fibersToKeep );

return