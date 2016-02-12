function [fe, fascicles_indices_in_new_fe] = feConnectomeReduceFibers(fe, fascicles_to_keep)
% Deletes a set of fibers from the signal and the M matrix.
%
%   [fe, fascicles_indices_in_new_fe] = feConnectomeDeleteFibers(fe,fascicles_to_keep)
%
% Inputs:
%   - fe, an fe structure, see feCreate.m, and v_lifeExample.m
%   - fascicles_to_keep, a list of indexes to the fibers to keep e.g., [1 10 100].
%   - reduceVoxels, (optional), if set to 1 it will remove all the voxels
%       from the fe structure where the fibers left inside the fe structure do
%       not have go through.
%
% Example:
%   fascicles_to_keep = 1:50;
%   feConnectomeReduceFibers(fe, fascicles_to_keep)
%
% Copyright (2015-2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% Find the indices of the fasciles kept in the new FE structure:
[fas_preserved, fas_i] = ismember(1:size(fe.life.M.Phi,3), fascicles_to_keep);
fascicles_indices_in_new_fe = fas_i(fas_preserved);

% Keeping the fibers we want in the indicator function Phi
fe.life.M.Phi = fe.life.M.Phi(:,:, (fascicles_to_keep) );

% Clear the fields that depend o the original fiber group. These are:
% (1) The fit of the model.
if isfield(fe.life,'fit') && ~isempty(fe.life.fit)
    fe.life.fit.weights         = fe.life.fit.weights(fascicles_to_keep);
    fe.life.fit.results.nParams = sum(fascicles_to_keep);
end
if isfield(fe.life,'voxfit') && ~isempty(fe.life.voxfit)
    fe.life.voxfit = [];
end
% The field 'fibers' containing some statistics obtained from the original
% fg
if isfield(fe.life,'fibers') && ~isempty(fe.life.fibers)
    fe.life.fibers.tensors = [];
    fe.life.fibers.total   = [];
    fe.life.fibers.unique  = [];
end

% Next, we will change the actual fiber group to reflect the virtual lesion.
if all(fascicles_to_keep==0)
    fe = feSet(fe,'fg img',fgCreate('all fibers were deleted',fe.fg.name,'img'));
else
    % feExtract requires indices to each fibers and does not accept logical
    % inputs. Here we make sure we are passing indices not a logical vector.
    if (length(unique(fascicles_to_keep)) == 2)
        if   (all(unique(fascicles_to_keep) == [0, 1]'))
            fascicles_to_keep = find(fascicles_to_keep);end
    end
    fg = feGet(fe,'fibers img');
    fg.fibers = fg.fibers(fascicles_to_keep);
    fe = feSet(fe,'fg img',fg);
end

return

