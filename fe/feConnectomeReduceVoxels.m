function [fe, indices_of_fas_in_new_fe] = feConnectomeReduceVoxels(fe,voxelsToKeep, fascicles_to_keep)
% Select the voxels to keep (extract) in a connectome matrix
%
%   [fe, indices_of_fas_in_new_fe] = feConnectomeReduceVoxels(fe,voxelToKeep)
%
%   voxelsToKeep:  Is a binary list of voxels we preserve.
%
% We expand the voxelsToKeep into a binary list of 0's and 1's that in
% which each is expanded by nBvecs.  The 1s are the rows of the M matrix we
% will keep.
%
%  Example:
%
% Copyright (2015-2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% Return only the mode and the signal for the voxels we want to keep
fe.life.M.Phi = fe.life.M.Phi(:,voxelsToKeep, :);

% Update the number of voxels inside the roi.
fe.roi.coords = fe.roi.coords(voxelsToKeep,:);

% Update the signal fields stored in the FE structure
fe.life.diffusion_signal_img  = fe.life.diffusion_signal_img(voxelsToKeep,:);
fe.life.diffusion_S0_img      = fe.life.diffusion_S0_img(voxelsToKeep);

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

% Remove the fascicles of the fascicle from the fe. We alos return the indices
% of all the fascicles we kept from the old FE structure into the new FE
% structure.
[fe, indices_of_fas_in_new_fe] = feConnectomeReduceFibers(fe, fascicles_to_keep );

return