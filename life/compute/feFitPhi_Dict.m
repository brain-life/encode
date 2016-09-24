function [Phi, Dict] = feFitPhi_Dict(varargin)
M = varargin{1};
w = varargin{2};
dSig = varargin{3};
fitMethod = varargin{4};
Niter = varargin{5};
bvecs = varargin{6};
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
        
        
        for iter=1:Niter
            B = ttv(M.Phi,w,3);
            [ind, val] = find(B);
            B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
            
            % Compute DB
            DB = M.DictSig*B;
            
            %% Update B(a,v) and Dict(theta,a)
            for a=1:nAtoms
                DB = DB - M.DictSig(:,a)*B(a,:);
                E = dSig - DB;
                pos = find(ind(:,1)==a);
                if ~isempty(pos)
                    cols = ind(pos,2);
                    
                    E = E(:,cols);
                    [da,ba] = right_nn_svds(E);
                    M.DictSig(:,a) = da;
                    B(a,cols) = ba;
                    
                    % uptdate DB
                    DB = DB + M.DictSig(:,a)*B(a,:);
                end
                
                %disp(['atom ',num2str(a),' of ',num2str(nAtoms)]);
            end
            disp(['Fit Dic and B, iter',num2str(iter),' error=',num2str(100*norm(dSig-M.DictSig*B,'fro')/norm(dSig,'fro')),' nnz Phi=',num2str(nnz(M.Phi))])
        end
        
        %% Compute Phi compatible with B(a,v)
        [sub, val] = find(M.Phi);
        ind = sub2ind(size(B), sub(:,1), sub(:,2));
        b = B(:);
        newval = w(sub(:,3)).*b(ind);
        newval = newval/sum(w.^2);
        
        Phi = sptensor(sub, newval, size(M.Phi));
        Dict = M.DictSig;
        
    otherwise
        error('Cannot fit tensor Phi using method: %s.\n',fitMethod);
end

end

function [da,ba] = right_nn_svds(E)
tol = 1e-8;
delta = Inf;

[u1,s1,v1] = svds(E,1);
% try positive sing vectors
da = u1;
ba = s1*v1;
ba(ba<0) = 0;
error = norm(E-da*ba','fro');
while delta > tol
    ba = da'*E;
    da = E*ba'/(ba*ba');
    errorn = norm(E-da*ba,'fro');
    delta = abs(errorn - error);
    error = errorn;
end
error_pos = errorn;
da_pos = da;
ba_pos = ba;

% try neg sing vectors
da = -u1;
ba = -s1*v1;
ba(ba<0) = 0;
error = norm(E-da*ba','fro');
while delta > tol
    ba = da'*E;
    da = E*ba'/(ba*ba');
    errorn = norm(E-da*ba,'fro');
    delta = abs(errorn - error);
    error = errorn;
end
error_neg = errorn;
da_neg = da;
ba_neg = ba;

if error_pos < error_neg
    da = da_pos;
    ba = ba_pos;
else
    da = da_neg;
    ba = ba_neg;
end

end
