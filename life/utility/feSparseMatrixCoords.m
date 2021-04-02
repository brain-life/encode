function [ out ] = feSparseMatrixCoords(M)
%[ mat ] = feSparseMatrixCoords(Phi)
%   Convert the ENCODE model object (Phi) and the dictionary signal (D) to 
%       a set of coordinates to be used in the 2d sparse matrix version of 
%       LiFE that incorporates streamline groupings during optimization.
%
% INPUTS:
%   M   - The model object from an initialized ENCODE object
%         example: M = feGet(fe, 'model';
%
% OUTPUT:
%   out - The conversion of the 3d sparse tensor and dictionary to the 
%         coordinates / values required for the 2d sparse matrix LiFE model
%
% Brent McPherson put this function together, but
% Andy Womack and Dan McDonald did most of the driving.
%
% Brent McPherson, Indiana University (c) 2021
%

% pull the number of streamlines
nfib = size(M.Phi, 3);

% Each slice (streamline) of Phi needs to be multiplied by the dictionary
% signal to create the 2d matrix for of each streamlines prediction. The
% streamlines can be stacked into a single 2d matrix by adding a column
% indicating which streamline the entries (subscripts, values) they
% correspond to.

disp([ 'Converting ' num2str(nfib) ' streamlines in Phi to sparse matrix coordinates...' ]);

% initialize the output
out = [];
%out = nan(size(M.Phi.vals,1), 4); % preallocate output to final size

% for every streamline
for fib = 1:nfib
    
    % pull the slice of the tensor
    slc = M.Phi(:,:,fib);
    
    % use the sptensor toolbox to multiply the slice by the dictionary for
    % making each streamlines prediction
    mm = ttm(slc, M.DictSig, 1);
    
    % create the index to identify the streamline in the new matrix
    idx = ones(size(mm.vals)) * ii;
    
    % for this streamline, store the values from the sparse product:
    % - x (dictionary) / y (voxel) subscripts into the sparse model matrix
    % - the index of the streamline being for this value
    % - the value stored at that entry
    iter = [ mm.subs, idx, mm.vals ];

    % append this streamline to the output
    out = [ out; iter ];
    %out(zzz, :) = iter; % need zzz to be indices of out for the streamline
    
end

%
% concerns / things to note
%
% mm.subs is type double from ttm - why?
% - this is dumb, these should be int - ought to fix
% - indices are integers, why store them as doubles?
% - BE CAREFUL WITH THE PRECISION WHEN SAVING THIS OUTPUT
% -- because they're saved as doubles, the integer indices will round if
%    the precision isn't high enough, resulting in errors.
%
% In theory, the sptensor toolbox should be able to transform the tensor
% to a matrix with a function distributed in the tool. However, that 
% doesn't appear to be the case (their documentation / notation is bad).
%
% Further, the tensor-matrix multiplication doesn't work for these types 
% (handling sparse objects sparsely). So many attempts to work on the 3d
% tensor directly resulted in comically large memory overflows. Which is
% weird, because they have a @sptensor/ttm function that ought to handle
% sparse tensors correctly. Maybe it worked in an older version of matlab?
%

end

