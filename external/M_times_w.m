function [y] = M_times_w(M,w)
% This function computes the matrix by vector product M*w
% using the factorized version of M which consists on only 
% two components: 
%       1) The Dictionary M.DictSig;
%       2) A sparse 3D array Phi with size [nAtoms,Nvoxels,nFibers]

% The following makes w'*Phi_{(3)} (multiplication in mode-1)
y = ttv(M.Phi,w,3); % This is memory efficient version of above

% The following makes ans \times_1 Dic
if nnz(y)
    y = ttm(y,M.DictSig,1);
else
    y = sptensor([size(M.DictSig,1),size(y,2)]); % In case empty tensor
end

%y = sptenmat(y,1); % sptensor -> sptenmat
y = tenmat(y,1); % sptensor -> tenmat
y = double(y); % tenmat -> sparse matrix
y = y(:); % vectorize

end

