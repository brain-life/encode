function [fe, removedFibers] = feConnectomeReduceFibers(fe, fibersToKeep)
% Deletes a set of fibers from the signal and the M matrix.
%
%   [fe, removedFibers] = feConnectomeDeleteFibers(fe,fibersToKeep)
%
% Inputs:
%   - fe, an fe structure, see feCreate.m, and v_lifeExample.m
%   - fibersToKeep, a list of indexes to the fibers to keep e.g., [1 10 100].
%   - reduceVoxels, (optional), if set to 1 it will remove all the voxels
%       from the fe structure where the fibers left inside the fe structure do
%       not have go through.
%
% Example:
%   fibersToKeep = 1:50;
%   feConnectomeReduceFibers(fe, fibersToKeep)
%
% Copyright (2015-2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% Keeping the fibers we want in the indicator function Phi
fe.life.M.Phi = fe.life.M.Phi(:,:, find(fibersToKeep) );

% Clear the fields that depend o the original fiber group. These are:
% (1) The fit of the model.
if isfield(fe.life,'fit') && ~isempty(fe.life.fit)
    fe.life.fit.weights         = fe.life.fit.weights(fibersToKeep);
    fe.life.fit.results.nParams = sum(fibersToKeep);
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
if all(fibersToKeep==0)
    fe = feSet(fe,'fg img',fgCreate('all fibers were deleted',fe.fg.name,'img'));
else
    % feExtract requires indices to each fibers and does not accept logical
    % inputs. Here we make sure we are passing indices not a logical vector.
    if (length(unique(fibersToKeep)) == 2)
        if   (all(unique(fibersToKeep) == [0, 1]'))
            fibersToKeep = find(fibersToKeep);end
    end
    fg = feGet(fe,'fibers img');
    fg.fibers = fg.fibers(fibersToKeep);
    fe = feSet(fe,'fg img',fg);
end

return

