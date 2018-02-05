function val = feGet(fe,param,varargin)
% Get function for fascicle evaluation structure
%
%   val = feGet(fe,param,varargin)
%
%
% INPUTS: 
%  Coords    - Nx3 set of coordinates in image space
%  voxelIndices - Vector of 1's and 0's, there is a one for each
%              location the the connectome coordinates for which there is
%              a match in the coords
%
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
%
%---------- List of arguments ----
% Name of the current fe structure.
% name = feGet(fe,'name')
%----------
% The type of objes (always, fascicle evaluation)
% type = feGet(fe,'type')
%----------
% Structure of parameters and results from LIFE analysis
% life = feGet(fe,'life')
%----------
% Return the connectome (fiber group) in image coordinates.
% fg = feGEt(fe,'fg img')
%----------
% return the connectome in acpc coordinates.
% fg  = feGet(fg,'fg acpc')
%----------
% Return the VOI comprised by the full connectome.
% roi = feGet(fe,'roi')
%----------
% Transforms between coords for fg, roi, and dwi data
% xform = feGet(fe,'xform')
%----------
% Load the diffusion weighted data.
% dwi = feGet(fe,'dwi')
%----------
% Return the path to the diffusion weighted data file.
% dwiFile = feGet(fe,'dwifile')   
%----------
% Load the repeated-measure diffusion weighted data.
% dwi = feGet(fe,'dwirepeat') 
%----------
% Load the repeated measure of the diffusion weighted data.
% dwiFile = feGet(fe,'dwirepeatfile')
%----------
% Directory where the LiFe structes are saved by defualt.
% sdir = feGet(fe,'savedir');
%----------
% Diffusion directions.
% val = feGet(fe,'bvecs');
%----------
% Indices to the diffusion directions in the DWi 4th Dimension.
% val = feGet(fe,'bvecs indices');
%----------
% Number of B0's
% val = feGet(fe,'n bvals');
%----------
% B0 Values.
% bval = feGet(fe,'bvals')
%----------
% Returns a nVoxels X nBvecs array of measured diffusion signal
% val = feGet(fe,'dsiinvox');
% val = feGet(fe,'dsiinvox',voxelsIndices);
% val = feGet(fe,'dsiinvox',coords);
%----------
% Returns a nVoxels X nBvecs array of demeaned diffusion signal
% val =feGet(fe,'dsiinvoxdemeaned');
% val =feGet(fe,'dsiinvoxdemeaned',voxelsIndices);
% val =feGet(fe,'dsiinvoxdemeaned',coords);
%----------
% Get the diffusion signal at 0 diffusion weighting (B0) for this voxel
% val = feGet(fe,'b0signalvoxel');
% val = feGet(fe,'b0signalvoxel',voxelIndex);
% val = feGet(fe,'b0signalvoxel',coords);
%----------
% Return a subset of fibers from the conenctome.
% fg = feGet(fe,'fiberssubset',fiberList);
%----------
% Return the coordinates of the VOI comprised in the connectome.
% Always in image space.
% coords = feGet(fe,'roi coords')
%----------
% Number of voxels in the connectome/VOI.
% nVoxels = feGet(fe,'n voxels')
%----------
% Indexes of actually used voxels.
% These can be different than the number of voxels in the VOI in case some
% voxels have no fibers in them.
%
% TO BE DEPRECATED.
% indxUsedVx = feGet(fe,'index of used voxels')
%----------
% Number of actually used voxels.
% This can be less than the number of voxels in the VOI in case some
% voxels have no fibers in them.
%
% TO BE DEPRECATED.
% nUsedVx = feGet(fe,'index of used voxels')
%----------
% Returns the first index to a voxel in the dSig vector.
% The dSig vector is 1:nBvecs*nVoxels.
% val = feGet(fe,'start rows')
%----------
% Return the rows corresponding to a set of voxels.
% rows = feGet(fe,'voxel rows',idxVoxel)
%----------
% Return the model (M matrix), or a subset of it.
% Mfiber = feGet(fe,'M fiber');
% Mfiber = feGet(fe,'model',voxelsIndices);
% Mfiber = feGet(fe,'model',coords);
%----------
% Return a subset of measurements from the fiber portion of M matrix
% The rows of M are specified.  Remember that the rows refer to a
% combination of voxel and measurement direction and potentially in the
% future a b-value.  Hence, the rows are not specific to a voxel.
% Mfiber = feGet(fe,'mfiber subset',rowsToKeep);
%----------
% Pairing of fibers and nodes in each voxel.
% v2fnp = feGet(fe,'v2fnp');
%----------
% Tensors computed for each node and fiber in a set of voxels.
% t = feGet(fe,'tensors')
% t = feGet(fe,'tensors',voxelsIndices);
% t = feGet(fe,'tensors',coords);
%----------
% Get all the tensors in a single voxel.
% val = feGet(fe,'voxtensors',voxelIndex)
% val = feGet(fe,'voxtensors',coord)
%----------
% Returns the model (M) with only a subset of columns (fibers).
% fiberList is a vector of indexes, e.g, [1 10 100]
% M = feGet(fe,'keepfibers',fiberList)
%----------
% Isotropic portion of M matrix. Returns the matrix for the full model or
% for a subset of voxels.
% Miso = feGet(fe,'M iso')
% Miso = feGet(fe,'M iso',voxelsIndices)
% Miso = feGet(fe,'M iso',coords)
%----------
% Total number of nodes in each voxel The voxel2FN pairs are row size
% equals number of nodes and column size is always 2. The entries of
% the first column are the fiber number. The entries in the second
% column are the node of the fiber inside the voxel.
% val = feGet(fe,'n nodes');
% val = feGet(fe,'n nodes',voxelsIndices);
% val = feGet(fe,'n nodes',coords);
%----------
% Return the total number of fibers for all the voxels or in a set of
% voxels.
% nFibers = feGet(fe,'uniquefnum');           % all the voxels
% nFibers = feGet(fe,'uniquefnum',[1 2 3 4]); % for some the voxels,
%                                             % specified by indexes
% nFibers = feGet(fe,'uniquefnum',coords);    % for some the voxels,
%                                             % specified by coordinates
%----------
% Return the indexes of the fibers for all the voxels or in a set of
% voxels.
% idxFibers = feGet(fe,'uniquef');           % all the voxels
% idxFibers = feGet(fe,'uniquef',[1 2 3 4]); % for some the voxels,
%                                            % specified by indexes
% idxFibers = feGet(fe,'uniquef',coords);    % for some the voxels,
%                                            % specified by coordinates
%----------
% Return the total number of fibers for all the voxels or in a set of
% voxels
% nFibers = feGet(fe,'totfnum');           % all the voxels
% nFibers = feGet(fe,'totfnum',[1 2 3 4]); % for some the voxels,
%                                          % specified by indexes
% nFibers = feGet(fe,'totfnum',coords);    % for some the voxels,
%                                          % specified by coordinates
%----------
% Return the indexes of the fibers for all the voxels or in a set of
% voxels
% idxFibers = feGet(fe,'totf');             % all the voxels
% idxFibers = feGet(fe,'totf',[1 2 3 4]);   % for some the voxels,
%                                           % specified by indexes
% idxFibers = feGet(fe,'totfibers',coords); % for some the voxels,
%                                           % specified by coordinates
%----------
% Return the number of fibers in the model.
% nFibers = feGet(fe,'n fibers');
%----------
% Weights of the isotropic voxel signals, this is the mean signal in
% each voxel.
% val = feGet(fe,'iso weights');
% val = feGet(fe,'iso weights'coords)
% val = feGet(fe,'iso weights',voxelIndices)
%----------
% Weights of the fiber component. For all the fibers or a subset of
% them.
% w = feGet(fe,'fiber weights')
% w = feGet(fe,'fiber weights',fiberIndices)
%----------
% Weights of the fiber component with a subset of the fibers' weights
% set to zero. This can be used to test the loss in RMSE for the
% connectome when a subset of fibers is removed, but the connectome is
% not fitted again.
% w = feGet(fe,'fiber weights test',fiberIndices)
%----------
% Predicted signal (demeaned) with a subset of fibers' weights set to 0.
% This can be used to test the loss in RMSE for the connectome when a
% subset of fibers is removed, but the connectome is not fitted again.
%
% pSig = feGet(fe,'pSig fiber test',fiberIndices);
% pSig = feGet(fe,'pSig fiber test',fiberIndices,voxelIndices);
% pSig = feGet(fe,'pSig fiber test',fiberIndices,coords);
%----------
% Measured signal in VOI, this is the raw signal. not demeaned
%
% dSig = feGet(fe,'dSig full')
%---------
% Measured signal in VOI, demeaned, this is the signal used for the
% fiber-portion of the M model.
%
% dSig = feGet(fe,'dsigdemeaned');
% dSig = feGet(fe,'dsigdemeaned',[1 10 100]);
% dSig = feGet(fe,'dsigdemeaned',coords);
%---------
% Get the demeaned signal for a subset of rows.
% Useful for cross-validation.
% dSig = feGet(fe,'dsigrowssubset',voxelsList);
%---------
% Predicted signal of fiber alone (demeaned).
% pSig = feGet(fe,'pSig fiber');
% pSig = feGet(fe,'pSig fiber',coords);
% pSig = feGet(fe,'pSig fiber',voxelIndices);
%---------
% Predicted measured signal from both fiber and iso
% pSig = feGet(fe,'pSig full')
% pSig = feGet(fe,'pSig full',coords)
% pSig = feGet(fe,'pSig full',voxelIndices)
%---------
% The woxels returned by a fit of LiFE by voxel/fiber
% w = feGet(fe,'fiberweightsvoxelwise')
%---------
% Predict the diffusion signal for the fiber component
% with the voxel-wise fit of LiFE
% pSig = feGet(fe,'psigfvoxelwise')
% pSig = feGet(fe,'psigfvoxelwise',coords)
% pSig = feGet(fe,'psigfvoxelwise',voxelIndices)
%---------
% Predict the diffusion signal for the fiber component
% with the voxel-wise fit of LiFE, return the an array of pSigXnVoxel
% pSig = feGet(fe,'psigfvoxelwisebyvoxel')
% pSig = feGet(fe,'psigfvoxelwisebyvoxel',coords)
% pSig = feGet(fe,'psigfvoxelwisebyvoxel',voxelIndices)
%---------
% The fiber and isotropic weights as a long vector
% w = feGet(fe,'fullweights')
%---------
% Return the global R2 (fraction of variance explained) of the full life
% model.
% R2 = feGet(fe,'total r2');
%---------
% Percent variance explained
% R2 = feGet(fe,'explained variance')
%---------
% Root mean squared error of the LiFE fit to the whole data
% rmse = feGet(fe,'rmse')
%---------
% Residual signal: (fiber prediction - measured_demeaned).
% res = feGet(fe,'res sig fiber')
%---------
% Residual signal: (full model prediction - measured signal).
% res = feGet(fe,'res sig full');
%---------
% Residual signal: (fiber model prediction - demeaned measured signal)
% with added mean signal. Res is returned as a vector.
% This is used to reconstruct an image (volume) to be used for the
% refinement process.
% res = feGet(fe,'fiber res sig with mean');
% res = feGet(fe,'fiber res sig with mean',coords);
% res = feGet(fe,'fiber res sig with mean',voxelIndices);
%---------
% Residual signal fiber model prediction - demeaned measured signal
% with added mean signal. Res is returned as an array (nBvecs x nVoxel).
% This is used to reconstruct an image (volume) to be used for the
% refinement process.
% res = feGet(fe,'fiber res sig with mean voxel');
% res = feGet(fe,'fiber res sig with mean voxel',coords);
% res = feGet(fe,'fiber res sig with mean voxel',voxelIndices);
%---------
% Predicted signal by the fiber model in a set of voxels.
% Vox subfield stores per voxel within the VOI. The pSig has size of the
% dwi. - It is stored as val = pSig(X,Y,Z,Theta)
%
% pSig = feGet(fe, 'pSig fiber by voxel');
% pSig = feGet(fe, 'pSig fiber by voxel',coords);
% pSig = feGet(fe, 'pSig fiber by voxel',voxelIndex);
%---------
% Return a column vector of the proportion of variance explained in
% each voxel.
% R2byVox = feGet(fe,'voxr2');
% R2byVox = feGet(fe,'voxr2',coords);
%---------
% Return a column vector of the proportion of variance explained in
% each voxel. (Normalized to the squared mean diffusion signal in each voxel)
% R2byVox = feGet(fe,'voxr2zero');
% R2byVox = feGet(fe,'voxr2zero',coords);
%---------
% Return the percent of varince explained in each voxel.
% R2byVox = feGet(fe,'var exp by voxel');
% R2byVox = feGet(fe,'var exp by voxel',coords);
%---------
% Demeaned diffusion signal in each voxel.
% dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel');
% dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',coords);
% dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',vxIndex);
%---------
% Full (measured) signal in each voxel.
% dSigByVoxel = feGet(fe,'dSig full by Voxel');
% dSigByVoxel = feGet(fe,'dSig full by Voxel',coords);
% dSigByVoxel = feGet(fe,'dSig full by Voxel',vxIndex);
%---------
% Predicted signal by the full model in a set of voxeles.
% pSigByVoxel = feGet(fe, 'pSig full by voxel');
% pSigByVoxel = feGet(fe, 'pSig full by voxel',coords);
% pSigByVoxel = feGet(fe, 'pSig full by voxel',voxelIndex);
%---------
% A volume of RMSE values.
% RMSE = feGet(fe,'vox rmse')
% RMSE = feGet(fe,'vox rmse',coords)
% RMSE = feGet(fe,'vox rmse',vxIndex)
%---------
% A volume of RMSE values with a subset of fibers' weights set to 0.
% RMSE = feGet(fe,'vox rmse',fiberIndices)
% RMSE = feGet(fe,'vox rmse',fiberIndices,coords)
% RMSE = feGet(fe,'vox rmse',fiberIndices,voxelIndex)
%---------
% Fibers' residual signal by voxel.
% res = feGet(fe,'res sig fiber vox')
% res = feGet(fe,'res sig fiber vox',coords)
% res = feGet(fe,'res sig fiber vox',vxIndex)
%---------
% Full (measured) residual signal by voxel.
% res = feGet(fe,'res sig full vox')
% res = feGet(fe,'res sig full vox',coords)
% res = feGet(fe,'res sig full vox',vxIndex)
%---------
% Fibers' residual signal by voxel with mean signal added (with added
% isotropic component).
% This is used to compute the residual signal for the refinement.
% res = feGet(fe,'fiber res sig with mean vox')
% res = feGet(fe,'fiber res sig with mean vox',coords)
% res = feGet(fe,'fiber res sig with mean vox',vxIndex)
%---------
% Residual signal full model prediction from a multi-voxel fit.
% res = feGet(fe,'res sig full voxfit');
% res = feGet(fe, 'res sig full voxfit',coords);
% res = feGet(fe, 'res sig full voxfit',voxelIndex);
%---------
% Fiber density statistics.
% Computes the fiber density (how many fibers in each voxel)
% We compute the following values:
% (1) The number of fibers in each voxel
% (2) The number of unique fibers with non-zero weights
% (3) The sum of the weights in each voxel
% (4) The mean of the weights in each voxel
% (5) The variance of the weigths in each voxel
% pSig = feGet(fe,'fiberdensity');
%---------
% Given an VOI or a set of indices to voxels returns the indices of the matching voxels inside the big
% volume, which ordinarily represents the full connectome.
% voxelIndices = feGet(fe,'voxelsindices',coords)
%                IMPORTANT: Size(coords) must be Nx3;
% voxelIndices = feGet(fe,'voxelsindices',voxelIndices)
%                IMPORTANT: Indices MUST be a column vector for this
%                code to work. Size(voxelIndices) = Nx1;
%---------
% Given an VOI finds the indices of the matching voxels inside the big
% volume, which ordinarily represents the full connectome.
%
% foundVoxels = feGet(fe,'find voxels',coords)
%
% coords is a Nx3 set of coordinates in image space
% foundVoxels is a vector of 1's and 0's, there is a one for each
% location the the connectome coordinates for which there is a match in
% the coords
%---------
% Given a set of VOI coords finds the row numbers of the Model matrix
% (or equivalently the dSig vector) that represent the data for this
% set of VOI coords.
% foundVoxels = feGet(fe,'coords 2 rows',coords)
% coords      - a Nx3 set of coordinates in image space
% foundVoxels - a binary vector that is 1 for each Model matrix row
%               that corresponds to at least one of the coords.
%---------
% Quaternian transformation from IMAGE space to ACPC space.
% xform = feGet(fe,'xform')
%---------
% Quaternian transformation from ACPC space to IMAGE space.
% xform = feGet(fe,'xform')
%---------
% Dimensions of the DW volume.
% dim = feGet(fe,'dims')
%---------
% Dimensions of the maps of parameters and results.
% dims = feGet(fe, 'mapsize')
%---------
% Path to the 3D Anatomy Volume.
% anatomyfile = feGet(fe, 't1 file')
%
% End of feGet.m parameters, 
   
val = [];

% Format the input parameters.
param = lower(strrep(param,' ',''));

% Start sorting the input and computing the output.
switch param
  case 'name'
    % Name of the current fe structure.
    %
    % name = feGet(fe,'name')
    val = fe.name;
    
  case 'type'
    % The type of objes (always, fascicle evaluation)
    %
    % type = feGet(fe,'type')
    val = fe.type; % Always fascicle evaluation
    
    % Set top level structure, not just single slot
  case 'life'
    % Structure of parameters and results from LIFE analysis
    %
    % life = feGet(fe,'life')
    val = fe.life; 
    
  case {'fgimg','fibersimg'}
    % Return the connectome (fiber group) in image coordinates.
    %
    % fg = feGet(fe,'fg img')
    val = fe.fg;  % Fiber group candidate fascicles, uses fgSet/Get
    
  case {'fgacpc','fibersacpc'}
    % return the connectome in acpc coordinates.
    %
    % fg  = feGet(fg,'fg acpc')
    xform = feGet(fe,'img2acpcxform');
    val   = dtiXformFiberCoords(feGet(fe,'fgimg'),xform,'acpc');
    
  case 'roi'
    % Return the VOI comprised by the full connectome.
    %
    % roi = feGet(fe,'roi')
    val = fe.roi;
    
  case 'roiacpc'
    % Return the VOI comprised by the full connectome.
    %
    % roi = feGet(fe,'roi acpc')    
    name   = sprintf('roi_%s',fe.fg.name); 
    randColor = rand(1,3);
    val = dtiNewRoi(name,randColor,fefgGet(feGet(fe,'fg acpc'),'unique image coords'));
    
  case 'xform'
    % Transforms between coords for fg, roi, and dwi data
    %
    % xform = feGet(fe,'xform')
    val = fe.life.xform;
    
  case 'dwi'
    %  Load the diffusion weighted data.
    %
    % dwi = feGet(fe,'dwi')
    val = dwiLoad(feGet(fe,'dwifile'));
  
  case 'dwirepeat'
    %  Load the diffusion weighted data.
    %
    % dwi = feGet(fe,'dwirepeat')
    val = dwiLoad(feGetRep(fe,'dwifile'));
    
  case 'dwifile'
    %  Load the diffusion weighted data.
    %
    % dwiFile = feGet(fe,'dwifile')
    val = fe.path.dwifile;
    
  case 'dwifilerep'
    %  Load the diffusion weighted data.
    %
    % dwiFile = feGet(fe,'dwifilerep')
    val = fe.path.dwifilerep;
 
  case 'savedir'
    %  Directory where the LiFe structes are saved by defualt.
    %
    % sdir = feGet(fe,'savedir');
    val = fe.path.savedir;
    
  case {'bvecs'}
    % Diffusion directions.
    %
    % val = feGet(fe,'bvecs');
    val = fe.life.bvecs;
    
  case {'bvecsindices'}
    % Indices to the diffusion directions in the DWi 4th Dimension.
    %
    % val = feGet(fe,'bvecs indices');
    val = fe.life.bvecsindices;
    
  case {'nbvecs','nbvals'}
    % Number of B0's
    %
    % val = feGet(fe,'n bvals');
    val = length(feGet(fe,'bvals'));
    
  case {'bvals'}
    % B0 Values.
    %
    % bval = feGet(fe,'bvals')
    val = fe.life.bvals;
    
  case {'diffusionsignalinvoxel','dsiinvox','dsigvox','dsigmeasuredvoxel'}
    % Returns a nVoxels X nBvecs array of measured diffusion signal
    %
    % val = feGet(fe,'dsiinvox');
    % val = feGet(fe,'dsiinvox',voxelsIndices);
    % val = feGet(fe,'dsiinvox',coords);
    val = fe.life.diffusion_signal_img(feGet(fe,'voxelsindices',varargin),:)';
    
  case {'diffusionsignalinvoxeldemeaned','dsiinvoxdemeaned'}
    % Returns a nVoxels X nBvecs array of demeaned diffusion signal
    %
    % val =feGet(fe,'dsiinvoxdemeaned');
    % val =feGet(fe,'dsiinvoxdemeaned',voxelsIndices);
    % val =feGet(fe,'dsiinvoxdemeaned',coords);
    nBvecs     = feGet(fe,'nBvecs');
    voxelIndices = feGet(fe,'voxelsindices',varargin);
    val = fe.life.diffusion_signal_img(voxelIndices,:) - repmat(mean(fe.life.diffusion_signal_img(voxelIndices,:), 2),1,nBvecs);
    keyboard
    % THis seems to be wrong
    
  case {'b0signalimage','b0vox'}
    % Get the diffusion signal at 0 diffusion weighting (B0) for this voxel
    %
    % val = feGet(fe,'b0signalvoxel');
    % val = feGet(fe,'b0signalvoxel',voxelIndex);
    % val = feGet(fe,'b0signalvoxel',coords);
    val = fe.life.diffusion_S0_img(feGet(fe,'voxelsindices',varargin), :);
    if length(size(val)) > 1
        val = mean(val,2);
    end
    
  case {'fiberssubset','fsub','subsetoffibers','fgsubset'}
    % Return a subset of fibers from the conenctome.
    %
    % fg = feGet(fe,'fiberssubset',fiberList);
    val = feGet(fe,'fg img');
    fiberList = varargin{1};
    fibers = cell(size(fiberList));
    for ff = 1:length(fiberList)
      fibers{ff} = val.fibers{fiberList(ff)};
      if isfield(val,'pathwayInfo')
        pathInfo(ff) = val.pathwayInfo(fiberList(ff));
      end
    end
    val.fibers = fibers;
    if isfield(val,'pathwayInfo')
      val.pathwayInfo = pathInfo;
    end
    val.name = sprintf('subsetFrom:%s',val.name);
    
  case {'roivoxels','roicoords'}
    % Return the coordinates of the VOI comprised in the connectome.
    % Always in image space.
    %
    % coords = feGet(fe,'roi coords')
    val = fe.roi.coords;
  
  case {'roicoordssubset'}
    % Return the coordinates of the VOI comprised in the connectome.
    % Always in image space.
    %
    % coords = feGet(fe,'roi coords')
    val = feGet(fe,'roi coords');
    val = val(varargin{1},:);
  
  case {'nroivoxels','nvoxels'}
    % Number of voxels in the connectome/VOI.
    %
    % nVoxels = feGet(fe,'n voxels')
    val = size(feGet(fe,'roivoxels'),1);
    
  case {'usedvoxels','indexesofusedvoxels'}
    % Indexes of actually used voxels.
    %
    % These can be different than the number of voxels in the VOI in case some
    % voxels have no fibers in them.
    %
    % TO BE DEPRECATED.
    %
    % indxUsedVx = feGet(fe,'index of used voxels')
    val = find(feGet(fe,'nnodes'));
    
  case {'nusedvoxels','numberofusedvoxels'}
    % Number of actually used voxels.
    %
    % This can be less than the number of voxels in the VOI in case some
    % voxels have no fibers in them.
    %
    % TO BE DEPRECATED.
    %
    % nUsedVx = feGet(fe,'index of used voxels')
    val = size(feGet(fe,'used voxels'),2);
    
  case {'startrow','first index into the dsig vector for a voxel'}
    % Returns the first index to a voxel in the dSig vector.
    %
    % The dSig vector is 1:nBvecs*nVoxels.
    %
    % val = feGet(fe,'start rows')
    val = (varargin{1}-1)*feGet(fe,'nBvecs') + 1;
    
  case {'voxelrows','all indexes into the dsig vector for a set of voxels'}
    % Return the rows corresponding to a set of voxels.
    %
    % rows = feGet(fe,'voxel rows',idxVoxel)
    n        = feGet(fe,'n voxels');
    nBvecs   = feGet(fe,'n bvecs');
    voxelsToKeep = zeros(n,1); voxelsToKeep(varargin{1}) = 1;
    val = logical(kron(voxelsToKeep(:),ones(nBvecs,1)));
    
  case {'modeltensor'}
    val = fe.life.modelTensor;
    
  case {'mfiber','m','model'}
    % Return the model (M matrix), or a subset of it.
    %
    % Mfiber = feGet(fe,'M fiber');
    % Mfiber = feGet(fe,'model',voxelsIndices);
    % Mfiber = feGet(fe,'model',coords);
    
    % idxVoxels is an integer list of the voxels to retain
    % If idxVoxels is not included, return the whole matrix
    % M contains only the rows for the specified voxels (all directions)
    
    val = fe.life.M;
        
  case {'voxel2fnpair','voxel2fibernodepair','v2fnp'}
    % Pairing of fibers and nodes in each voxel.
    %
    % v2fnp = feGet(fe,'v2fnp');
    if isempty(fe.life.voxel2FNpair)
      fprintf('[%s] fe.life.voxel2FNpair is empty, can be computed as: \nfe = feSet(fe,''v2fnp'',feGet(fe,''fg img''),''v2fn'',feGet(fe,''roi coords'')',mfilename);
      return
    end
    val = fe.life.voxel2FNpair;
    
  case {'tensors','fibers tensors'}
    % Tensors computed for each node and fiber in a set of voxels.
    %
    % t = feGet(fe,'tensors')
    % t = feGet(fe,'tensors',voxelsIndices);
    % t = feGet(fe,'tensors',coords);
    if isempty(varargin)
      val = fe.life.fibers.tensors;
    else
      vv = feGet(fe,'voxelsindices',varargin);
      val = cell(size(vv));
      for ff = 1:length(vv)
        val{ff} = fe.life.fibers.tensors{vv(ff)};
      end
    end
    
  case {'voxeltensors','voxtensors','voxq'}
    % Get all the tensors in a signle voxel
    %
    % val = feGet(fe,'voxtensors',voxelIndex)
    % val = feGet(fe,'voxtensors',coord)
    
    % Index for the voxel
    vv = feGet(fe,'voxelsindices',varargin);
    
    % The indexes of the voxeles used, to build LiFE.
    usedVoxels = feGet(fe,'usedVoxels');
    voxIndex   = usedVoxels(vv);
    nNodes     = feGet(fe,'nNodes');
    
    % Get the tensors for each node in each fiber going through this voxel:
    val = zeros(nNodes(voxIndex), 9); % Allocate space for all the tensors (9 is for the 3 x 3 tensor components)
    for ii = 1:nNodes(voxIndex)           % Get the tensors
      val(ii,:) = fe.life.fibers.tensors{fe.life.voxel2FNpair{voxIndex}(ii,1)} ...
        (fe.life.voxel2FNpair{voxIndex}(ii,2),:);
    end
    
  case {'nnodes','numofnodes'}
    % Total number of nodes in each voxel The voxel2FN pairs are row size
    % equals number of nodes and column size is always 2. The entries of
    % the first column are the fiber number. The entries in the second
    % column are the node of the fiber inside the voxel.
    %
    % val = feGet(fe,'n nodes');
    % val = feGet(fe,'n nodes',voxelsIndices);
    % val = feGet(fe,'n nodes',coords);
    if ~isempty(fe.life.voxel2FNpair)
        [val, ~] = cellfun(@size,fe.life.voxel2FNpair);
        val      = val(feGet(fe,'voxelsindices',varargin));
    end

  case {'numberofuniquefibersbyvoxel','uniquefnum'}
    % Return the total number of fibers for all the voxels or in a set of
    % voxels
    %
    % nFibers = feGet(fe,'uniquefnum');           % all the voxels
    % nFibers = feGet(fe,'uniquefnum',[1 2 3 4]); % for some the voxels,
    %                                             % specified by indexes
    % nFibers = feGet(fe,'uniquefnum',coords);    % for some the voxels,
    %                                             % specified by coordinates
    val = fe.life.fibers.unique.num(feGet(fe,'voxelsindices',varargin));
    
  case {'indextouniquefibersbyvoxel','uniquef'}
    % Return the indexes of the fibers for all the voxels or in a set of
    % voxels.
    %
    % idxFibers = feGet(fe,'uniquef');           % all the voxels
    % idxFibers = feGet(fe,'uniquef',[1 2 3 4]); % for some the voxels,
    %                                            % specified by indexes
    % idxFibers = feGet(fe,'uniquef',coords);    % for some the voxels,
    %                                            % specified by coordinates
    vxIndex = feGet(fe,'voxels indices',varargin);
    val = cell(length(vxIndex),1);
    for ii = 1:length(vxIndex)
      val{ii} = fe.life.fibers.unique.index{vxIndex(ii)};
    end
    
  case {'numberoftotalfibersbyvoxels','totfnum'}
    % Return the total number of fibers for all the voxels or in a set of
    % voxels
    %
    % nFibers = feGet(fe,'totfnum');           % all the voxels
    % nFibers = feGet(fe,'totfnum',[1 2 3 4]); % for some the voxels,
    %                                          % specified by indexes
    % nFibers = feGet(fe,'totfnum',coords);    % for some the voxels,
    %                                          % specified by coordinates
    val = fe.life.fibers.total.num(feGet(fe,'voxels indices',varargin));
    
  case {'indexoftotalfibersbuvoxel','totf'}
    % Return the indexes of the fibers for all the voxels or in a set of
    % voxels
    %
    % idxFibers = feGet(fe,'totf');             % all the voxels
    % idxFibers = feGet(fe,'totf',[1 2 3 4]);   % for some the voxels,
    %                                           % specified by indexes
    % idxFibers = feGet(fe,'totfibers',coords); % for some the voxels,
    %                                           % specified by coordinates
    vxIndex = feGet(fe,'voxels indices',varargin);
    val = cell(length(vxIndex),1);
    for ii = 1:length(vxIndex)
      val{ii} = fe.life.fibers.total.index{vxIndex(ii)};
    end
    
  case {'nfibers'}
    % Return the number of fibers in the model.
    %
    %val = size(fe.life.fibers,2);
    val = size(fe.fg.fibers,1);
    
  case {'natoms'}
    % Return the number of atoms in the dictionary
    val = size(fe.life.M.DictSig,2);  

  case {'orient'}
    % Return the number of atoms in the dictionary
    val = fe.life.M.orient;     
    
  case {'isoweights','weightsiso','meanvoxelsignal'}
    % Weights of the isotropic voxel signals, this is the mean signal in
    % each voxel.
    %
    % val = feGet(fe,'iso weights');
    % val = feGet(fe,'iso weights'coords)
    % val = feGet(fe,'iso weights',voxelIndices)
    val = feGet(fe,'Miso') \ feGet(fe,'dsig full')';
    val = val(feGet(fe,'voxelsindices',varargin));
     
  case {'voxisoweights','voxweightsiso','voxmeanvoxelsignal'}
    % Weights of the isotropic voxel signals, this is the mean signal in
    % each voxel.
    % Organized nBvecs x nVoxels
    %
    % val = feGet(fe,'iso weights');
    % val = feGet(fe,'iso weights'coords)
    % val = feGet(fe,'iso weights',voxelIndices)
    val = feGet(fe,'Miso') \ feGet(fe,'dsig full')';
    val = repmat(val,1,feGet(fe,'nbvecs'))';

    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(:,feGet(fe,'voxelsindices',varargin));
    end

  case {'dsigiso','isodsig'}
    % mean signalin each voxel, returned for each diffusion direction:
    % size(nBvecsxnVoxels, 1).
    %
    % val = feGet(fe,'iso disg');
    % val = feGet(fe,'iso disg',coords)
    % val = feGet(fe,'iso dsig',voxelIndices)
    val = feGet(fe,'iso weights');
    val = repmat(val,1,feGet(fe,'nbvecs'))';
    val = val(:);
    
  case {'fiberweights'}
    % Weights of the fiber component. For all the fibers or a subset of
    % them.
    %
    % w = feGet(fe,'fiber weights')
    % w = feGet(fe,'fiber weights',fiberIndices)

    % If the model was not fit yet, fit it, install the fit and then return
    % the weights.
    if ~isfield(fe.life,'fit')
      fprintf('[%s] fe.life.fit is empty, can be computed as: \nfeSet(fe,''fit'',feFitModel(feGet(fe,''Mfiber''), feGet(fe,''dsigdemeaned''),''sgdnn''))',mfilename);
      return
    end
    
    % Get the weights
    if ~isempty(varargin) % subselect the weights for the requested fibers
      val = fe.life.fit.weights(varargin{1});  
    else % Return all of them
      val = fe.life.fit.weights;
    end
    
  case {'fiberweightstest'}
    % Weights of the fiber component with a subset of the fibers' weigth
    % set to zero. This can be used to test the loss in RMSE for the
    % connectome when a subset of fibers is removed, but the connectome is
    % not fitted again.
    %
    % w = feGet(fe,'fiber weights test',fiberIndices)

    % If the model was not fit yet, fit it, install the fit and then return
    % the weights.
    if ~isfield(fe.life,'fit')
      fprintf('[%s] fe.life.fit is empty, can be computed as: \nfeSet(fe,''fit'',feFitModel(feGet(fe,''Mfiber''), feGet(fe,''dsigdemeaned''),''sgdnn''))',mfilename);
      return
    end
    
    % Get the weights
    if ~isempty(varargin) % subselect the weights for the requested fibers
      val = fe.life.fit.weights;  
      val(varargin{1}) = 0;
    else % A set of fiber weights is required as input    
      error('[%s] Indices to a subset of fiber-weights must be passed in w = feGet(fe,''fiber weights test'',fiberIndices))',mfilename);
    end
    
    case {'fiberweightstestvoxelwise'}
    % Weights of the fiber component with a subset of the fibers' weigth
    % set to zero. This can be used to test the loss in RMSE for the
    % connectome when a subset of fibers is removed, but the connectome is
    % not fitted again.
    %
    % w = feGet(fe,'fiber weights test voxel wise',fiberIndices)

    % If the model was not fit yet, fit it, install the fit and then return
    % the weights.
    if ~isfield(fe.life,'voxfit')
      fprintf('[%s] fe.life.fit is empty, can be computed as: \nfeSet(fe,''fit'',feFitModel(feGet(fe,''Mfiber''), feGet(fe,''dsigdemeaned''),''sgdnn''))',mfilename);
      return
    end
    
    % Get the weights
    if ~isempty(varargin) % subselect the weights for the requested fibers
      val = fe.life.voxfit.weights;  
      val(:,varargin{1}) = 0;
    else % A set of fiber weights is required as input    
      error('[%s] Indices to a subset of fiber-weights must be passed in\nw = feGet(fe,''fiber weights test voxel wise'',fiberIndices))',mfilename);
    end
    
  case {'dsigmeasured','dsigfull'}
    % Measured signal in VOI, this is the raw signal. not demeaned
    %
    % dSig = feGet(fe,'dSig full')
    val = fe.life.diffusion_signal_img';
    % Return a subset of voxels
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val              = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'diffusionsignaldemeaned','dsigdemeaned'}
    % Measured signal in VOI, demeaned, this is the signal used for the
    % fiber-portion of the M model.
    %
    % dSig = feGet(fe,'dsigdemeaned');
    % dSig = feGet(fe,'dsigdemeaned',[1 10 100]);
    % dSig = feGet(fe,'dsigdemeaned',coords);
    nVoxels = feGet(fe,'nVoxels');
    nBvecs  = feGet(fe,'nBvecs');
    dSig = reshape(fe.life.diffusion_signal_img',[1,nVoxels*nBvecs]);
    val     = (dSig - reshape(repmat( ...
      mean(reshape(dSig, nBvecs, nVoxels),1),...
      nBvecs,1), size(dSig)))';
    % Return a subset of voxels
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end

  case {'dsigrowssubset','diffusionsignaldemeanedinsubsetofrows'}
    % Get the demeaned signal for a subset of rows.
    % Useful for cross-validation.
    %
    % dSig = feGet(fe,'dsigrowssubset',voxelsList);
    dSig = feGet(fe,'dsigdemeaned');
    for vv = 1:length(voxelsList)
      val(feGet(fe,'voxel rows',vv)) = dSig(feGet(fe,'voxel rows',voxelList(vv)));
    end
    
  case {'psigfiber', 'fiberpsig','fiberpredicted','dsigpredictedfiber'}
    % Predicted signal of fiber alone (demeaned).
    %
    % pSig = feGet(fefeGet(fe,'fiber weights','pSig fiber');
    % pSig = feGet(fe,'pSig fiber',coords);
    % pSig = feGet(fe,'pSig fiber',voxelIndices);
    %val = feGet(fe,'Mfiber')*feGet(fe,'fiber weights');
    %val = M_times_w(feGet(fe,'Mfiber'),feGet(fe,'fiber weights'));
    nTheta  = feGet(fe,'nbvecs');
    nVoxels = feGet(fe,'nvoxels');  
    val = M_times_w(fe.life.M.Phi.subs(:,1),fe.life.M.Phi.subs(:,2),fe.life.M.Phi.subs(:,3),fe.life.M.Phi.vals,fe.life.M.DictSig,feGet(fe,'fiber weights'),nTheta,nVoxels);
    
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'uniquefibersindicesinroi'}
    % Find the unique fibers indices in the FE roi.
    % Thisis necessary sometimes after changing the number of voxels in an
    % FE structure, as it is performed by feConnectomeReduceVoxels.m
    %
    % FibInRoi = feGet(fe,'uniquefibersinroi');
    
    % Get all the unique fibers in each voxel
    uniquefvx = fefgGet(feGet(fe,'fibers img'), ...
                        'uniquefibersinvox',     ...
                        feGet(fe,'roi coords'));
    val = [];
    for ivx = 1:length(uniquefvx)
        val = [val; uniquefvx{ivx}];
    end
    val = unique(val);
    
  case {'weightsinroi'}
    % Find the weights of the fibers in a specified volume.
    %
    % We compute the following values:   
    % (1) The number of unique fibers in the volume
    % (2) The weights for the fibers in the roi
    %
    % wFibInRoi = feGet(fe,'weightsinroi');
    
    % Get all the unique fibers in each voxel
    FibInRoi = feGet(fe,'uniquefibersindicesinroi');
   
    % extract the fber weights obtained in a LiFE fit
    w   = feGet(fe,'fiber weights');
    
    % Return only the ones for the fibers in this roi
    val = w(unique( FibInRoi ));
        
  case {'nnzw','numberofnonzerofibers','connectomedensity'}
    % Connectome density (number of non-zero weighted fascicles).
    %
    % conDensity = feGet(fe,'nnz');
                          
    % extract the fber weights obtained in a LiFE fit
    fw = feGet(fe,'fiber weights');
    val = sum(fw > 0);
    
  case {'indnzw'}
    %
    % Return the indices of the nonzero weights
    %
    fw = feGet(fe, 'fiber weights');
    val = find(fw > 0 & ~isnan(fw));
    
  case {'fiberdensity'}
    % Fiber density statistics.
    %
    % Computes the fiber density (how many fibers in each voxel)
    % 
    % We compute the following values:   
    % (1) The number of fibers in each voxel
    % (2) The number of unique fibers with non-zero weights
    % (3) The sum of the weights in each voxel
    % (4) The mean of the weights in each voxel
    % (5) The variance of the weigths in each voxel
    %
    % pSig = feGet(fe,'fiberdensity');
    
    % Get the unique fibers in each voxel
    uniquefvx = fefgGet(feGet(fe,'fibers img'), ...
                        'uniquefibersinvox',     ...
                        feGet(fe,'roi coords'));
                      
    % extract the fber weights obtained in a LiFE fit
    w = feGet(fe,'fiber weights');
    
    % Compute the fiber density in three wasy:
    % (1) The number of fibers in each voxel
    % (2) The number of unique fibers with non-zero weights
    % (3) The sum of the weights in each voxel
    % (4) The mean of the weights in each voxel
    % (5) The variance of the weigths in each voxel
    val = nan(length(uniquefvx),5);
    for ivx = 1:length(uniquefvx)
        
      % Number of fibers in each voxel
      val(ivx,1) = length(uniquefvx{ivx});
      
      % Number of fibers in ech voxel with non-zero weight
      val(ivx,2) = length(uniquefvx{ivx}(w(uniquefvx{ivx}) > 0));
  
      % Sum of fiber weights in each voxel
      val(ivx,3) = sum(w(uniquefvx{ivx}));
         
      % Mean of fiber weights in each voxel
      val(ivx,4) = nanmedian(w(uniquefvx{ivx}));
      
      % Var of fibers in each voxel
      val(ivx,5) = nanvar(w(uniquefvx{ivx}));
    end
    
  case {'fullweights'}
    % The fiber and isotropic weights as a long vector
    %
    % w = feGet(fe,'fullweights')
    val = [feGet(fe,'fiber weights'); feGet(fe,'iso weights')];
    
  case {'psigfvoxelwisebyvoxel'}
    % Predict the diffusion signal for the fiber component 
    % with the voxel-wise fit of LiFE, return the an array of pSigXnVoxel
    %
    % pSig = feGet(fe,'psigfvoxelwisebyvoxel')
    % pSig = feGet(fe,'psigfvoxelwisebyvoxel',coords)
    % pSig = feGet(fe,'psigfvoxelwisebyvoxel',voxelIndices)
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'psigfvoxelwise'), nBvecs, nVoxels);

  case {'totalr2'}
    % Return the global R2 (fraction of variance explained) of the full life
    % model.
    %
    % R2 = feGet(fe,'total r2');

    %     measured  = feGet(fe,'dsigdemeaned');
    %     predicted = feGet(fe,'fiber p sig');
    %     val = (1 - (sum((measured - predicted).^2 ) ./ ...
    %                 sum((measured - mean(measured)).^2) ));
    val = (1 - (sum((feGet(fe,'diffusion signal demeaned') - ...
      feGet(fe,'fiber predicted')).^2 ) ./ ...
      sum((feGet(fe,'diffusion signal demeaned') - ...
      mean(feGet(fe,'diffusion signal demeaned'))).^2) ));
   
  case {'totpercentvarianceexplained','totpve'}
    % Percent variance explained
    %
    % R2 = feGet(fe,'explained variance')
    val = 100 * feGet(fe,'total r2');
    
  case {'rmsetotal','totalrmse'}
    % Root mean squared error of the LiFE fit to the whole data
    %
    % rmse = feGet(fe,'rmse')
    val = sqrt(mean((feGet(fe,'diffusion signal demeaned') - ...
      feGet(fe,'pSig fiber')).^2));
  
  case {'totalrmsevoxelwise'}
    % Root mean squared error of the LiFE fit to the whole data from a
    % vocel-wise fit
    %
    % rmse = feGet(fe,'rmse')
    val = sqrt(mean((feGet(fe,'diffusion signal demeaned') - ...
      feGet(fe,'pSig f voxelwise')').^2));
  
  case {'ressigfiber'}
    % Residual signal fiber prediction - measured_demeaned.
    %
    % res = feGet(fe,'res sig fiber')
    val = (feGet(fe,'dsigdemeaned')    - feGet(fe,'psig fiber'));
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val              = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'ressigfull'}
    % Residual signal full model prediction - measured signal.
    %
    % res = feGet(fe,'res sig full');
    % res = feGet(fe, 'res sig full',coords);
    % res = feGet(fe, 'res sig full',voxelIndex);
    val = (feGet(fe,'dsig full')    - feGet(fe,'psig full')');
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'ressigfullvoxelwise'}
    % Residual signal full model prediction from a multi-voxel fit.
    %
    % res = feGet(fe,'res sig full voxfit');
    % res = feGet(fe, 'res sig full voxfit',coords);
    % res = feGet(fe, 'res sig full voxfit',voxelIndex);
    %tic,val      = feGet(fe,'dsig demeaned') - feGet(fe,'psigfvoxelwise') + feGet(fe,'dsigiso');toc
    val      = feGet(fe,'dsig full') - feGet(fe,'psigfvoxelwise')';
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'fiberressigwithmean'}
    % Residual signal fiber model prediction - demeaned measured signal
    % with added mean signal. VECTOR FORM.
    %
    % This is used to reconstruct an image (volume) to be used for the
    % refinement process.
    %
    % res = feGet(fe,'fiber res sig with mean');
    
    % predicted = feGet(fe,'psig fiber');
    % measured  = feGet(fe,'dsigdemeaned');
    % val = (measured_full - predicted demeaned);
    val = (feGet(fe,'dsig full')' - feGet(fe,'psig fiber')); 
    
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'predictedfibersignalvoxel','psigfvox'}
    % Predicted signal by the fiber model in a set of voxels.
    %
    % Vox subfield stores per voxel within the VOI
    % The pSig has size of the dwi.
    % It is stored as val = pSig(X,Y,Z,Theta)
    %
    % pSig = feGet(fe, 'pSig fiber by voxel');
    % pSig = feGet(fe, 'pSig fiber by voxel',coords);
    % pSig = feGet(fe, 'pSig fiber by voxel',voxelIndex);
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'pSig fiber'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin));
    
  case {'predictedfullsignalvoxel','psigfullvox'}
    % Predicted signal by the full model in a set of voxeles.
    %
    % Vox subfield stores per voxel within the VOI
    % The pSig has size of the dwi.
    % It is stored as val = pSig(X,Y,Z,Theta)
    %
    % pSig = feGet(fe, 'pSig full by voxel');
    % pSig = feGet(fe, 'pSig full by voxel',coords);
    % pSig = feGet(fe, 'pSig full by voxel',voxelIndex);
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'pSig full'), nBvecs, nVoxels);
      
  case {'voxelr2','r2vox','voxr2','r2byvoxel'}
    % Return a column vector of the proportion of variance explained in
    % each voxel.
    %
    % R2byVox = feGet(fe,'voxr2');
    % R2byVox = feGet(fe,'voxr2',coords);
    measured  = feGet(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    nBvecs    = feGet(fe,'nBvecs');
    val       = (1 - (sum((measured - predicted).^2 ) ./ ...
      sum((measured - repmat(mean(measured), nBvecs,1)).^2 ) ));
    val       = val(feGet(fe,'return voxel indices',varargin));
  
  case {'voxelr2zero','r2voxzero','voxr2zero','r2byvoxelzero'}
    % Return a column vector of the proportion of variance explained in
    % each voxel. (Normalized to the squared diffusion signal in each voxel)
    %
    % R2byVox = feGet(fe,'voxr2zero');
    % R2byVox = feGet(fe,'voxr2zero',coords);
    measured  = feGet(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    val       = (1 - (sum((measured - predicted).^2 ) ./ sum(measured.^2)));
    val       = val(feGet(fe,'return voxel indices',varargin));
   
  case {'voxelvarianceexplained','varexpvox','voxvarexp','varexpbyvoxel'}
    % Return the percent of varince explained in each voxel.
    %
    % R2byVox = feGet(fe,'var exp by voxel');
    % R2byVox = feGet(fe,'var exp by voxel',coords);
    val = 100.*feGet(fe,'voxr2');
    val = val(feGet(fe,'return voxel indices',varargin));
    
  case {'voxelr2voxelwise','r2voxvoxelwise','voxr2voxelwise'}
    % Return a column vector of the proportion of variance explained in
    % each voxel.
    %
    % R2byVox = feGet(fe,'voxr2voxelwise');
    % R2byVox = feGet(fe,'voxr2voxewise',coords);
    measured  = feGet(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'psig f voxel wise by voxel');
    nBvecs    = feGet(fe,'nBvecs');
    val       = (1 - (sum((measured - predicted).^2 ) ./ ...
      sum((measured - repmat(mean(measured), nBvecs,1)).^2 ) ));
    val       = val(feGet(fe,'return voxel indices',varargin));
  
  case {'voxelr2zerovoxelwise','voxr2zerovoxelwise'}
    % Return a column vector of the proportion of variance explained in
    % each voxel. (Normalized to the squared diffusion signal in each voxel)
    %
    % R2byVox = feGet(fe,'voxr2zerovoxelwise');
    % R2byVox = feGet(fe,'voxr2zerovoxelwise',coords);
    measured  = feGet(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'psig f voxel wise by voxel');
    val       = (1 - (sum((measured - predicted).^2 ) ./ sum(measured.^2)));
    val       = val(feGet(fe,'return voxel indices',varargin));
   
  case {'voxelvarianceexplainedvoxelwise','varexpvoxvoxelwise'}
    % Return the percent of varince explained in each voxel.
    %
    % R2byVox = feGet(fe,'var exp by voxel');
    % R2byVox = feGet(fe,'var exp by voxel',coords);
    val = 100.*feGet(fe,'voxr2voxelwise');
    val = val(feGet(fe,'return voxel indices',varargin));
    
  case {'dsigdemeanedbyvoxel','dsigdemeanedvox','Y'}
    % Demeaned diffusion signal in each voxel
    %
    % dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel');
    % dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',coords);
    % dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',vxIndex);
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'dsigdemeaned'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin)); 
    
  case {'dsigfullbyvoxel','dsigfullvox','voxdsigfull'}
    % Full (measured) signal in each voxel
    %
    % dSigByVoxel = feGet(fe,'dSig full by Voxel');
    % dSigByVoxel = feGet(fe,'dSig full by Voxel',coords);
    % dSigByVoxel = feGet(fe,'dSig full by Voxel',vxIndex);
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'dSig full'), nBvecs, nVoxels);
    val    = val(:,feGet(fe,'return voxel indices',varargin))';
    
  case {'voxelrmse','voxrmse'}
    % A volume of RMSE values
    %
    % RMSE = feGet(fe,'vox rmse')
    % RMSE = feGet(fe,'vox rmse',coords)
    % RMSE = feGet(fe,'vox rmse',vxIndex)
    measured  = feGet(fe,'dsigdemeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    val       = sqrt(mean((measured - predicted).^2,1));
    val       = val(feGet(fe,'voxelsindices',varargin));
   
  case {'voxrmses0norm'}
      % A volume of RMSE normalized by the S0 value in each voxel.
      rmse = feGet(fe,'vox rmse');
      s0   = feGet(fe,'b0signalimage')';
      % Some voxels can have a S0=0. We replace the S0 in these voxles with
      % a NaN. 
      idx = (s0 == 0);
      rmse(idx) = nan(size(find(idx)));
      
      val = rmse./s0;
        
  case {'voxelrmsevoxelwise','voxrmsevoxelwise'}
    % A volume of RMSE values optained with the voxel-wise (voxelwise) fit.
    %
    % RMSE = feGet(fe,'vox rmse voxelwise')
    % RMSE = feGet(fe,'vox rmse voxelwise',coords)
    % RMSE = feGet(fe,'vox rmse voxelwise',vxIndex)
    measured  = feGet(fe,'dsigdemeaned by voxel');
    predicted = feGet(fe,'pSig f voxel wise by voxel');
   
    val       = sqrt(mean((measured - predicted).^2,1));
    val       = val(feGet(fe,'voxelsindices',varargin));

  case {'voxelrmsetest','voxrmsetest'}
    % A volume of RMSE values with a subset of fibers' weights set to 0.
    %
    % RMSE = feGet(fe,'vox rmse',fiberIndices)
    % RMSE = feGet(fe,'vox rmse',fiberIndices,coords)
    % RMSE = feGet(fe,'vox rmse',fiberIndices,voxelIndex)
    measured  = feGet(fe,'dsigdemeaned by voxel');
    % Reshape the predicted signal by voxles
    predicted = reshape(feGet(fe,'pSig fiber test',varargin{1}), feGet(fe,'nBvecs'), feGet(fe,'nVoxels'));   
    val       = sqrt(mean((measured - predicted).^2,1));
    
    if length(varargin) == 2
      val       = val(feGet(fe,'voxelsindices',varargin));
    end
    
  case {'residualsignalfibervoxel','resfibervox'}
    % Fibers' residual signal by voxel
    %
    % res = feGet(fe,'res sig fiber vox')
    % res = feGet(fe,'res sig fiber vox',coords)
    % res = feGet(fe,'res sig fiber vox',vxIndex)
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'res sig fiber'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
    
  case {'residualsignalfullvoxel','resfullvox'}
    % Full (measured) residual signal by voxel
    %
    % res = feGet(fe,'res sig full vox')
    % res = feGet(fe,'res sig full vox',coords)
    % res = feGet(fe,'res sig full vox',vxIndex)
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'res sig full'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
      
  case {'fiberressigwithmeanvoxel'}
    % Residual signal fiber model prediction - demeaned measured signal
    % with added mean signal. VOXEL FORM (nBvecs x nVoxel).
    %
    % This is used to reconstruct an image (volume) to be used for the
    % refinement process.
    %
    % res = feGet(fe,'fiber res sig with mean');
    
    % predicted = feGet(fe,'psig fiber');
    % measured  = feGet(fe,'dsigdemeaned');
    % val = (measured     - predicted) + (measured_full - measured);
    val     = feGet(fe,'fiberressigwithmean');
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(val, nBvecs, nVoxels);
    
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val              = val(voxelRowsToKeep,:);
      val = val(:, feGet(fe,'return voxel indices',varargin));
    end
    
  case {'fiberressigwithmeanvox','resfibermeanvox'}
    % Fibers' residual signal by voxel with mean signal added (with added
    % isotropic component).
    %
    % This is used to compute the residual signal for the refinement.
    %
    % res = feGet(fe,'fiber res sig with mean vox')
    % res = feGet(fe,'fiber res sig with mean vox',coords)
    % res = feGet(fe,'fiber res sig with mean vox',vxIndex)
    nBvecs  = feGet(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGet(fe,'fiber res sig with mean'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
    
  case {'voxelsindices','returnvoxelindices','voxelindexes''returnvoxelindexes'}
    % Given an VOI or a set of indices to voxels returns the indices of the matching voxels inside the big
    % volume, which ordinarily represents the full connectome.
    %
    % voxelIndices = feGet(fe,'voxelsindices',coords)
    % voxelIndices = feGet(fe,'voxelsindices',voxelIndices)
    %                IMPORTANT: Indices MUST be a column vector for this
    %                code to work.
    %
    % coords is a Nx3 set of coordinates in image space
    % voxelIndices is a vector of 1's and 0's, there is a one for each
    % location the the connectome coordinates for which there is a match in
    % the coords
    %
    varargin = varargin{1};
    if ( ~isempty(varargin) )
      if ( size(varargin{1},2) == 3 )
        % If a set of coordinates were passed in we return the indices.
        val = logical(feGet(fe,'find voxels',varargin{1}));
      else
        % If indices were passed in, we sort and return them.
        val = sort(varargin{1});
      end
    else
      % If no coordinates were passed at all, we return all the indices.
      val = true(feGet(fe,'n voxels'),1);
    end
    
  case {'findvoxels','findvoxelindexes','findvoxelsinconnectome','findvoxelsinconnectomeroi', ...
      'voxel2index','coords2index'}
    % Given an VOI finds the indices of the matching voxels inside the big
    % volume, which ordinarily represents the full connectome.
    %
    % foundVoxels = feGet(fe,'find voxels',coords)
    %
    % coords is a Nx3 set of coordinates in image space
    % foundVoxels is a vector of 1's and 0's, there is a one for each
    % location the the connectome coordinates for which there is a match in
    % the coords
    
    % Compute the values for a subset of voxels
    if isempty(varargin)
      % If no coordinates were passed at all, we return all the indices.
      val = ones(feGet(fe,'n voxels'),1);
    else
      % The stored coordinates are at a resolution of the image, typically
      % 2mm isotropic.  The VOI should also be in image resolution.
      % ---- Franco:
      % This is how I had it. I think it is worng:
      % val = ismember(feGet(fe,'roi coords'), varargin{1}, 'rows'); % This is slow
      % This is hwo I think it should be:
      [~,val] = ismember(varargin{1}, feGet(fe,'roi coords'),'rows'); % This is slow
    end
    
  case {'voxelcoords2voxelrows','coords2rows'}
    % Given a set of VOI coords finds the row numbers of the Model matrix
    % (or equivalently the dSig vector) that represent the data for this
    % set of VOI coords.
    %
    % foundVoxels = feGet(fe,'coords 2 rows',coords)
    %
    % coords      - a Nx3 set of coordinates in image space
    % foundVoxels - a binary vector that is 1 for each Model matrix row
    %               that corresponds to at least one of the coords.
    val = feGet(fe,'voxel rows',find(feGet(fe,'find voxels',varargin{1})));
    
    % ------ Spatial coordinate transforms for voxels and fg to coordinate frames
  case {'xformimg2acpc','xform2acpc','img2acpc','img2acpcxform'}
    % Quaternian transformation from IMAGE space to ACPC space.
    %
    % xform = feGet(fe,'xform')
    val = fe.life.xform.img2acpc;
  case {'xformacpc2img','xform2img','acpc2img','acpc2imgxform'}
    % Quaternian transformation from ACPC space to IMAGE space.
    %
    % xform = feGet(fe,'xform')
    val = fe.life.xform.acpc2img;
  case {'volumesize','dims','dim','imagedim'}
    % Dimensions of the DW volume.
    %
    % dim = feGet(fe,'dims')
    val = fe.life.imagedim;
  case {'mapsize'}
    % Dimensions of the maps of parameters and results.
    %
    % dims = feGet(fe, 'mapsize')
    val = fe.life.imagedim(1:3);
    
  case {'anatomyfile'}
    % Path to the 3D Anatomy Volume.
    % anatomyfile = feGet(fe, 'anatomy file')
    val = fe.path.anatomy;
  case 'dtfile'
    % Diffusion weighted file used for testing results
    % dtfile = feGet(fe,'dtfile')
     val = fe.path.dtfile;  
  case 'phi'
    % sparse core tensor
     val = fe.life.M.Phi;
  case {'dictionary', 'd'}
     val = fe.life.M.DictSig;
        
  case 'lists'
    % atoms, voxels and Nelem lists (to be revised for future elimination)
     val.atoms = fe.life.fibers.atoms_list;       
     val.voxels = fe.life.fibers.voxels_list;       
     val.Nelem = fe.life.fibers.Nelem;  
  case 's0_img'
     val = fe.life.diffusion_S0_img;
  case 'nelem'
     val = fe.life.fibers.Nelem;
  case 'nphi'
     val = fe.life.M.Nphi;     
  case 'ntheta'
     val = fe.life.M.Ntheta;    
  case 'pathneighborhood'
      %
      % Find the path neighborhood of a set of fascicles. In this case the
      % set of fascicles is intepreted as a 'white matter tract' and the
      % fascicles sharing voxels with the tract are called path
      % neighborhood and their indices returned.
      %
      % INPUT: indices of fascicles composing a tract.
      % OUTPUT: indices of the fascicles sharing the same voxels.
      
      % provides the indices of fibers that touch the same voxels where a
      % provided tract exists
      [inds, ~] = find(fe.life.M.Phi(:,:,varargin{1})); % find nnz entries of subtensor      
           % inds has a list of (i,j,k) positions of nnz entries. Since we are interested in
      % locating voxels we need to look at the second column (j).
      voxel_ind = unique(inds(:,2)); 
      
      % To find which other fibers crosses these voxels, we need to look at
      % the the subtensor that corresponds to those voxels
      % See following lines
      [inds, ~] = find(fe.life.M.Phi(:,voxel_ind,:)); % find indices for nnz in the subtensor defined by voxel_ind
      val       = unique(inds(:,3)); % Fibers are the 3rd dimension in the subtensor
      val       = setdiff(val,varargin{1});
      
      % Find nnz weights indices and filter the obtained path neighborhood
      %w       = feGet(fe,'fiber weights');
      %ind_nnz = find(w);
      if isfield(fe.life.fit,'weights')
        w       = feGet(fe,'fiber weights');
        ind_nnz = find(w);
      else
        ind_nnz = unique(fe.life.M.Phi.subs(:,3));
      end
      val     = intersect(ind_nnz,val);
      
    case 'indnnz'
        val = unique(fe.life.M.Phi.subs(:,3));
        
    case 'voxindfromfibers'
        disp('Serching roi from fibers ...');
        [inds, ~] = find(fe.life.M.Phi(:,:,varargin{1})); % find nnz entries of subtensor
        val = unique(inds(:,2));
        
    case 'coordsfromfibers'
        disp('Serching roi from fibers ...');
        [inds, ~] = find(fe.life.M.Phi(:,:,varargin{1})); % find nnz entries of subtensor
        voxel_ind = unique(inds(:,2));
        val = feGet(fe,'roicoords');
        val = val(voxel_ind,:);
    case 'relativeerror'
        dSig = feGet(fe,'dsigdemeaned');
        dSig_pred = feGet(fe,'psigfiber');
        val = norm(dSig - dSig_pred)/norm(dSig);
        
    case 'predfull'
        [nAtoms]    = feGet(fe,'natoms');
        [nFibers]   = feGet(fe,'nfibers');
        [nVoxels]   = feGet(fe,'nvoxels');
        [nTheta]    = feGet(fe,'nbvals');
        D = fe.life.M.Dict;
        Phi = fe.life.M.Phi;
        
        if isfield(fe.life.fit, 'weights')
            B = ttv(Phi,fe.life.fit.weights,3);
        else
            B = ttv(Phi,ones(nFibers,1),3);
        end
        [ind, val] = find(B);
        B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        s0 = fe.life.s0;
        
        val = ones(nTheta,1)*s0' + D*B; 
        
    case 'predtract'
        [nAtoms]    = feGet(fe,'natoms');
        [nFibers]   = feGet(fe,'nfibers');
        [nVoxels]   = feGet(fe,'nvoxels');
        [nTheta]    = feGet(fe,'nbvals');
        D = fe.life.M.Dict;
        Phi = fe.life.M.Phi;
        
        w = zeros(nFibers,1);
        if isfield(fe.life.fit, 'weights')
            w(varargin{1}) = fe.life.fit.weights(varargin{1});
        else
            w(varargin{1})=1;
        end
        B = ttv(Phi,w,3);
        [ind, val] = find(B);
        B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        
        %s0 = fe.life.s0;
        %val = ones(nTheta,1)*s0' + D*B; % with iso   
        val = D*B; % without iso 
        
     case 'predtractwithiso'
        [nAtoms]    = feGet(fe,'natoms');
        [nFibers]   = feGet(fe,'nfibers');
        [nVoxels]   = feGet(fe,'nvoxels');
        [nTheta]    = feGet(fe,'nbvals');
        D = fe.life.M.Dict;
        Phi = fe.life.M.Phi;
        
        w = zeros(nFibers,1);
        w(varargin{1})=1;
        B = ttv(Phi,w,3);
        [ind, val] = find(B);
        B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        
        s0 = fe.life.s0;
        val = ones(nTheta,1)*s0' + D*B; % with iso  

     case 'predlesionwithiso'
        [nAtoms]    = feGet(fe,'natoms');
        [nFibers]   = feGet(fe,'nfibers');
        [nVoxels]   = feGet(fe,'nvoxels');
        [nTheta]    = feGet(fe,'nbvals');
        D = fe.life.M.Dict;
        Phi = fe.life.M.Phi;
        
        w = ones(nFibers,1);
        w(varargin{1})=0;
        B = ttv(Phi,w,3);
        [ind, val] = find(B);
        B = sparse(ind(:,1),ind(:,2),val,nAtoms,nVoxels);
        
        s0 = fe.life.s0;
        val = ones(nTheta,1)*s0' + D*B; % with iso          
        
  otherwise
    help('feGet')
    fprintf('[feGet] Unknown parameter << %s >>...\n',param);
    keyboard
end

end % END MAIN FUNCTION
