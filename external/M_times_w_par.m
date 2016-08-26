function [y] = M_times_w_par(M,w)
% This function computes the matrix by vector product M*w
% using the factorized version of M which consists on only 
% two components: 
%       1) The Dictionary M.DictSig;
%       2) A sparse 3D array Phi with size [nAtoms,Nvoxels,nFibers]

% The following makes w'*Phi_{(3)} (multiplication in mode-1)


Phi{1} = M.Phi(:,1:2500,:);
Phi{2} = M.Phi(:,2501:5000,:);
Phi{3} = M.Phi(:,5001:7500,:);
Phi{4} = M.Phi(:,7501:end,:);

tic
spmd
    y = ttv(Phi{labindex},w,3);    
end
disp(['Time with spmd:',num2str(toc)]);

subs = [];
vals = [];
for block = 1:4
    a = y{block};
    newsubs = a.subs;
    newsubs(:,2) = newsubs(:,2) + 2500*(block-1); 
    subs = vertcat(subs,newsubs);
    vals = vertcat(vals,a.vals);
end

y = sptensor(subs,vals,[32221,10694]);



tic
spmd(1)
    yn = ttv(M.Phi,w,3); % This is memory efficient version of above
end
disp(['Time WITHOUT spmd:',num2str(toc)]);

% 
% tic
% yn2 = ttv(Phi{1},w,3); % This is memory efficient version of above
% toc
% 
% tic
% yn2 = ttv(Phi{2},w,3); % This is memory efficient version of above
% toc
% 
% tic
% yn2 = ttv(Phi{3},w,3); % This is memory efficient version of above
% toc
% 
% tic
% yn2 = ttv(Phi{4},w,3); % This is memory efficient version of above
% toc


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

