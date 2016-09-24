function [Phi] = feFitPhi(varargin)
M = varargin{1};
w = varargin{2};
dSig = varargin{3};
fitMethod = varargin{4};
Niter = varargin{5};
%preconditioner = varargin{5};


% feFitPhi() function that given a fix vector of weights w optimize the tensor Phi
%   
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%load variables.mat

switch fitMethod
    case {'bbnnls'}
        [nFibers] = size(M.Phi,3); %feGet(fe,'nfibers');
        [nAtoms] = size(M.DictSig,2); %feGet(fe,'natoms');
        [nTheta] = size(M.DictSig,1);
        [nVoxels] = size(M.Phi,2); %feGet(fe,'nvoxels');
        
        
        dSig = reshape(dSig,[nTheta,nVoxels]);
        B = ttv(M.Phi,w,3);
        [ind, val] = find(B);
        B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        
        
        %% Update B(a,v)
        %options = optimset('Display','final','TolFun',1e-8,'TolX',1e-4);  
        options = optimset('lsqnonneg');
        for v=1:nVoxels
            bv = B(:,v);
            if nnz(bv)
                [ind, j, val] = find(bv);
                
                %result = bbnnls_orig(M.DictSig(:,ind), dSig(:,v), zeros(length(ind),1), opt);
                %x = result.x;
                [x, resnorm, residual, exitflag, output] = lsqnonneg(M.DictSig(:,ind), dSig(:,v),options);
                %x(x==0) = eps;
                bvn = sparse(ind,j,x,nAtoms,1);
                
                difbv = norm(bv-bvn)/norm(bv);
                B(ind,v) = x;
                error_pre = norm(dSig(:,v) - M.DictSig*bv)/norm(dSig(:,v));
                error = norm(dSig(:,v) - M.DictSig*bvn)/norm(dSig(:,v));
                %error = sqrt(resnorm)/norm(dSig(:,v));
                
                %disp(['Voxel ',num2str(v),' Error Pre=',num2str(100*error_pre),'%'' Error Post=',num2str(100*error),'%',' iterations=',num2str(output.iterations),' dif bv=',num2str(difbv)])
            end
        end
        
        %% Compute Phi compatible with B(a,v)
        [sub, val] = find(M.Phi);
        ind = sub2ind(size(B), sub(:,1), sub(:,2));
        b = B(:);
        newval = w(sub(:,3)).*b(ind);
        newval = newval/sum(w.^2);
        
        Phi = sptensor(sub, newval, size(M.Phi));

        
    otherwise
        error('Cannot fit tensor Phi using method: %s.\n',fitMethod);
end

end
