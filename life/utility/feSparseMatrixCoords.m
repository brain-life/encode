function [ out ] = feSparseMatrixCoords(M)
%[ mat ] = feSparseMatrixCoords(Phi)
%   Convert the ENCODE model object (Phi) and the dictionary signal (D) to 
%       a set of coordinates to be used in the 2d sparse matrix version of 
%       LiFE that incorporates streamline groupings during optimization.
%
% INPUTS:
%   M   - The model object from an initialized ENCODE object
%         example: M = feGet(fe, 'model');
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

fprintf('Converting %d streamlines in Phi to sparse matrix coordinates...\n', nfib);

% initialize the output
%out = [];
out = cell(nfib, 1);

% for every streamline
for fib = 1:nfib
    
    % print out every 1000th streamline as an update
    if mod(fib, 1000) == 0
        fprintf('Running streamline %d\n', fib);
    end
    
    % pull the slice of the model tensor for the streamline
    slc = M.Phi(:,:,fib);
    
    % use the sptensor toolbox to multiply the slice by the dictionary for
    % making each streamlines prediction
    mm = ttm(slc, M.DictSig, 1);
    
    % create the index to identify the streamline in the new matrix
    idx = ones(size(mm.vals)) * fib;
    
    % for this streamline, store the values from the sparse product:
    % - x (dictionary) / y (voxel) subscripts into the sparse model matrix
    % - the index of the streamline for these values
    % - the value stored at each entry
    iter = [ mm.subs, idx, mm.vals ];
    
    % append this streamline to the output matrix
    %out = [ out; iter ];
    out{fib} = iter;
    
end

fprintf('Stacking the streamline data into sparse matrix coordinates...\n');

% merge the output
out = cat(1, out{:});

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
% (handling sparse objects sparsely). So attempts to work on the 3d tensor 
% directly resulted in comically large memory overflows. Which is weird, 
% because they have a @sptensor/ttm function that ought to handle sparse 
% tensors correctly. Maybe it worked in an older version of matlab?
%

end

