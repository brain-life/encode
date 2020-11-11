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
fprintf(['\n[%s] Computing demeaned and non-demeaned difussivities dictionaries in a (',num2str(Nphi),'x',num2str(Ntheta),')-grid', ' ...'],mfilename); 
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

Dict = zeros(nBvecs,Norient); % Initialize Signal Dictionary matrix
DictSig = zeros(nBvecs,Norient); % Initialize Signal Dictionary matrix
%DictTensors = zeros(9,Norient); % Initialize Tensors Dictionary matrix

D = diag(fe.life.modelTensor); % diagonal matix with diffusivities
% THIS WILL GO IN THE LOOP BELOW

%[ shells, ~, shellindex ] = unique(round(bvals)); % round to get around noise in shells - FRAGILE FIX

ubv = feGet(fe, 'nshells');
ubi = feGet(fe, 'shellindex');

% Compute each dictionary column for a different kernel orientation
for j=1:Norient
    [Rot,~, ~] = svd(orient(:,j)); % Compute the eigen vectors of the kernel orientation
    
    Q = Rot*D*Rot'; % THIS GETS PULLED BY SHELL
    %DictTensors(:,j) = Q(:);
    Dict(:,j) = exp(- bvals .* diag(bvecs*Q*bvecs')); % Compute the signal contribution of a fiber in the kernel orientation divided S0
    % THIS CAN BE MODIFIED - HARD-CODED TENSOR - COULD BE DKI?
    
    % for every shell - THIS INCORPORATES EVERYTHING BUT ROT
    for k=1:size(ubv, 1) 
        si = ubi == ubv(k); % find the indices for the shell
        DictSig(si,j) = Dict(si,j) - median(Dict(si,j)); % demedianed signal by shell (used to demean)
    end
    
end

fe = feSet(fe,'dictionary parameters',{Nphi,Ntheta,orient,Dict,DictSig});

end