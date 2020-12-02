function fe = feBuildDictionaries(fe,Nphi,Ntheta)
% - Ntheta: number of discretization steps in Theta (azimuth)
% - Nphi: number of discretization steps in Phi (pi/2 - elevetion)
% - orient: (3 x Norient) containing in its columns the unit
% vectors covering half of the 3D unit sphere. These are the orientations
% that will be used to compute the dictionary of kernels.
% Norient = (Nphi-2)(Ntheta-2) + 1. The first column in orient is a vector
% pointing out to the zenith (spin-up) the rest of columns are vectors
% covering the half sphere.
%
%  Copyright (2020) Indiana University
%
%  Franco Pestilli frakkopesto@gmail.com and 
%  Cesar F. Caiafa ccaiafa@gmail.com

tic
fprintf(['\n[%s] Computing demeaned and non-demeaned diffusivities dictionaries in a (',num2str(Nphi),'x',num2str(Ntheta),')-grid', ' ...'],mfilename); 
fprintf('took: %2.3fs.\n',toc)

% Compute orientation vectors
Norient = (Ntheta-1)*Nphi + 1;
orient = zeros(3,Norient);
deltaTheta = pi/Ntheta;
deltaPhi = pi/Nphi;

sinTheta = sin(deltaTheta:deltaTheta:pi-deltaTheta);
cosTheta = cos(deltaTheta:deltaTheta:pi-deltaTheta);
sinPhi = sin(0:deltaPhi:pi-deltaPhi);
cosPhi = cos(0:deltaPhi:pi-deltaPhi);

orient(:,1) = [0;0;1];

orient(1,2:end) = kron(cosPhi,sinTheta);
orient(2,2:end) = kron(sinPhi,sinTheta);
orient(3,2:end) = repmat(cosTheta,[1,Nphi]);

% Compute Dictionary of Demeaned Signals for Canonical Diffusivities
nBvecs       = feGet(fe,'nBvecs');
bvecs        = feGet(fe,'bvecs');                      % bvecs
bvals        = feGet(fe,'bvals');                      % bvals

%bvecs(:,1) = -bvecs(:,1); % better observed x predicted correlation
%bvals = bvals * 10; % does bad

Dict = zeros(nBvecs,Norient); % Initialize Signal Dictionary matrix
DictSig = zeros(nBvecs,Norient); % Initialize Signal Dictionary matrix
%DictTensors = zeros(9,Norient); % Initialize Tensors Dictionary matrix

% catch Kurtosis estimates for debugging
%KDict = zeros(nBvecs,Norient); 
%KDictSig = zeros(nBvecs,Norient); 
akc = zeros(nBvecs,1);

% pull the shell information
ubv = feGet(fe, 'nshells');
ubi = feGet(fe, 'shellindex');
ubl = unique(ubi);

% pull diffusion tensor
dt = feGet(fe, 'model tensor');

% pull kurtosis fit - needs feGet
kt = fe.life.modelKurtosis;

% build apparent kurtosis coefficient - akc
for i=1:ubv
    
    % pull tensor parameters for hard indexing of equations
    sdt = dt(i,:);
    
    % find the indices for the shell
    si = ubi == ubl(i); 
        
    % mean diffusivity of tensor from forward model
    md = mean(sdt);
        
    % apparent diffusion coefficient - matches dipy
    adc = bvecs(si,1) .* bvecs(si,1) * sdt(1) + ...
          2 * bvecs(si,1) .* bvecs(si,2) * sdt(2) + ...
          bvecs(si,2) .* bvecs(si,2) .* sdt(3) + ...
          2 * bvecs(si,1) .* bvecs(si,3) * 0 + ... % only store primary eigenvalues in
          2 * bvecs(si,2) .* bvecs(si,3) * 0 + ... % LiFE prediction. These are 0.
          bvecs(si,3) .* bvecs(si,3) * 0;
    
    % apparent diffusion variance - matches dipy
    adv = ...
          bvecs(si,1) .* bvecs(si,1) .* bvecs(si,1) .* bvecs(si,1) .* kt(1) + ...       % xxxx
          bvecs(si,2) .* bvecs(si,2) .* bvecs(si,2) .* bvecs(si,2) .* kt(2) + ...       % yyyy
          bvecs(si,3) .* bvecs(si,3) .* bvecs(si,3) .* bvecs(si,3) .* kt(3) + ...       % zzzz
          4 * bvecs(si,1) .* bvecs(si,1) .* bvecs(si,1) .* bvecs(si,2) .* kt(4) + ...   % xxxy
          4 * bvecs(si,1) .* bvecs(si,1) .* bvecs(si,1) .* bvecs(si,3) .* kt(5) + ...   % xxxz
          4 * bvecs(si,1) .* bvecs(si,2) .* bvecs(si,2) .* bvecs(si,2) .* kt(6) + ...   % xyyy
          4 * bvecs(si,1) .* bvecs(si,3) .* bvecs(si,3) .* bvecs(si,3) .* kt(7) + ...   % xzzz
          4 * bvecs(si,2) .* bvecs(si,2) .* bvecs(si,2) .* bvecs(si,3) .* kt(8) + ...   % yyyz
          4 * bvecs(si,2) .* bvecs(si,3) .* bvecs(si,3) .* bvecs(si,3) .* kt(9) + ...   % yzzz
          6 * bvecs(si,1) .* bvecs(si,1) .* bvecs(si,2) .* bvecs(si,2) .* kt(10) + ...  % xxyy
          6 * bvecs(si,1) .* bvecs(si,1) .* bvecs(si,3) .* bvecs(si,3) .* kt(11) + ...  % xxzz
          6 * bvecs(si,2) .* bvecs(si,2) .* bvecs(si,3) .* bvecs(si,3) .* kt(12) + ...  % yyzz
          12 * bvecs(si,1) .* bvecs(si,1) .* bvecs(si,2) .* bvecs(si,3) .* kt(13) + ... % xxyz
          12 * bvecs(si,1) .* bvecs(si,2) .* bvecs(si,2) .* bvecs(si,3) .* kt(14) + ... % xyyz
          12 * bvecs(si,1) .* bvecs(si,2) .* bvecs(si,3) .* bvecs(si,3) .* kt(15);      % xyzz
    
    % zero out noisy (bad) adc estimates
    adc(adc < 0) = 0;
    
    % estimate apparent kutosis coefficient - matches dipy
    akc(si) = adv .* ((md ./ adc).^2);
    
    % zero out noisy (bad) akc estimates
    akc(akc < -3/7) = -3/7;
        
end

% Compute each dictionary column for a different kernel orientation
for j=1:Norient
    
    % Compute the eigen vectors of the kernel orientation
    [Rot,~, ~] = svd(orient(:,j));
    
    % for every shell
    for k=1:ubv
        
        % find the indices for the shell
        si = ubi == ubl(k); 
        
        % create diagonal matix with diffusivities for current shell
        % this assumes tensor fits for shell are entered in the order this will parse them in
        D = diag(dt(k,:));
        
        % estimate Q for tensor values in shell
        Q = Rot*D*Rot';
        
        % Compute the signal contribution of a fiber in the kernel orientation divided S0
        Dict(si,j)  = exp(-bvals(si) .* diag(bvecs(si,:)*Q*bvecs(si,:)')); 
        %Dict(si,j) = exp(-bvals(si) .* diag(bvecs(si,:)*Q*bvecs(si,:)') + (-bvals(si).^2 .* diag(bvecs(si,:)*Q*bvecs(si,:)').^2 .* akc(si))/6);
        
        % demeaned signal by shell
        DictSig(si,j) = Dict(si,j) - mean(Dict(si,j)); 
        %KDictSig(si,j) = KDict(si,j) - mean(KDict(si,j));
        
    end
    
end

fe = feSet(fe,'dictionary parameters',{Nphi,Ntheta,orient,Dict,DictSig});

% hard add kurtosis dictionaries / akc for debugging
%fe.life.M.KDict = KDict;
%fe.life.M.KDictSig = KDictSig;
fe.life.M.akc = akc; 

end