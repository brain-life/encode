function [ w ] = Mtransp_times_b(atoms,voxels,fibers,values,D,Y,nFibers)
% function [ w ] = Mtransp_times_b( M, b )
% [nFibers] = size(M.Phi,1);
% [nTheta]  = size(M.DictSig,1);
% [nAtoms] = size(M.Phi,2);
% [Nvoxels] = size(M.Phi,3);
% 
% w = M.DictSig'*reshape(b,[nTheta,Nvoxels]); % This is still a little bit 
% % memory expensive because w is a (Nd x Nvoxels) matrix (Nd = # atoms)
% % See comments below for a very memory efficient implementation but very
% % slow without using a compiled version (MEX).
% M.Phi = reshape(M.Phi,[nFibers,nAtoms*Nvoxels]);
% w = double(ttv(M.Phi,w(:),2));
% end
%
%  
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
% The following is a memory efficient version that should be compiled in
% order to provide fast results
w = zeros(nFibers,1);
for k = 1:length(values)
    %k
    w(fibers(k)) = w(fibers(k)) + (D(:,atoms(k)))'*Y(:,voxels(k))*values(k);
end

% 
% for f=1:Nf
%     for l=1:Nelem(f)
%         w(f) = w(f) + (D(:,atoms(f,l)))'*Y(:,voxels(f,l))*S0(voxels(f,l));
%     end
% end


% [nFibers] = size(M.Phi,1);
% [nTheta]  = size(M.DictSig,1);
% [Nvoxels] = size(M.Phi,3);
% 
% b = reshape(b,[nTheta,Nvoxels]); % matrix reshape
% 
% w = zeros(nFibers,1); % output
% 
% for f=1:nFibers
%     f
%     A = sptenmat(M.Phi(f,:,:),1); % ->sptenmat
%     aux = sum(b(:,A.subs(:,2)).*M.DictSig(:,A.subs(:,1)),1);
%     w(f) = sum(aux'.*A.vals(:));
% end
% 
% 
% end

