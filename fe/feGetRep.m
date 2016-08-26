function val = feGetRep(fe,param,varargin)
% Get function for fascicle evaluation structure, to use for a repeated
% measure calculation. 
%
% This function is similar to feGet but uses the repeted dataset to compute
% the DWI signal and all the calculations that depend on the DWI signal
% will be different.
%
% Importantly, this function will use calls to feGet every time we compute
% values that depend on the LiFE model or fit.
%
% This function allows for computing measures of:
% (1) data reliability and 
% (2) model versus data reliability
%
%   val = feGetRep(fe,param,varargin)
%
% INPUTS: 
%  Coords    - Nx3 set of coordinates in image space
%  voxelIndices - Vector of 1's and 0's, there is a one for each
%              location the the connectome coordinates for which there is
%              a match in the coords
%
%
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
%---------- List of arguments ----
% Name of the current fe structure.
% name = feGetRep(fe,'name')
%----------
% The type of objes (always, fascicle evaluation)
% type = feGetRep(fe,'type')
%----------
% Load the repeated-measure diffusion weighted data.
% dwi = feGetRep(fe,'dwirepeat') 
%----------
% Load the repeated measure of the diffusion weighted data.
% dwiFile = feGetRep(fe,'dwirepeatfile')
%----------
% Directory where the LiFe structes are saved by defualt.
% sdir = feGetRep(fe,'savedir');
%----------
% Diffusion directions.
% val = feGetRep(fe,'bvecs');
%----------
% Indices to the diffusion directions in the DWi 4th Dimension.
% val = feGetRep(fe,'bvecs indices');
%----------
% Number of B0's
% val = feGetRep(fe,'n bvals');
%----------
% B0 Values.
% bval = feGetRep(fe,'bvals')
%----------
% Returns a nVoxels X nBvecs array of measured diffusion signal
% val = feGetRep(fe,'dsiinvox');
% val = feGetRep(fe,'dsiinvox',voxelsIndices);
% val = feGetRep(fe,'dsiinvox',coords);
%----------
% Returns a nVoxels X nBvecs array of demeaned diffusion signal
% val =feGetRep(fe,'dsiinvoxdemeaned');
% val =feGetRep(fe,'dsiinvoxdemeaned',voxelsIndices);
% val =feGetRep(fe,'dsiinvoxdemeaned',coords);
%----------
% Get the diffusion signal at 0 diffusion weighting (B0) for this voxel
% val = feGetRep(fe,'b0signalimage');
% val = feGetRep(fe,'b0signalimage',voxelIndex);
% val = feGetRep(fe,'b0signalimage',coords);
%----------
% Weights of the isotropic voxel signals, this is the mean signal in
% each voxel.
% val = feGetRep(fe,'iso weights');
% val = feGetRep(fe,'iso weights'coords)
% val = feGetRep(fe,'iso weights',voxelIndices)
%----------
% Measured signal in VOI, this is the raw signal. not demeaned
%
% dSig = feGetRep(fe,'dSig full')
%---------
% Measured signal in VOI, demeaned, this is the signal used for the
% fiber-portion of the M model.
%
% dSig = feGetRep(fe,'dsigdemeaned');
% dSig = feGetRep(fe,'dsigdemeaned',[1 10 100]);
% dSig = feGetRep(fe,'dsigdemeaned',coords);
%---------
% Get the demeaned signal for a subset of rows.
% Useful for cross-validation.
% dSig = feGetRep(fe,'dsigrowssubset',voxelsList);
%---------
% Return the global R2 (fraction of variance explained) of the full life
% model.
% R2 = feGetRep(fe,'total r2');
%---------
% Percent variance explained
% R2 = feGetRep(fe,'explained variance')
%---------
% Root mean squared error of the LiFE fit to the whole data
% rmse = feGetRep(fe,'rmse')
%---------
% Residual signal: (fiber prediction - measured_demeaned).
% res = feGetRep(fe,'res sig fiber')
%---------
% Residual signal: (full model prediction - measured signal).
% res = feGetRep(fe,'res sig full');
%---------
% Residual signal: (fiber model prediction - demeaned measured signal)
% with added mean signal. Res is returned as a vector.
% This is used to reconstruct an image (volume) to be used for the
% refinement process.
% res = feGetRep(fe,'fiber res sig with mean');
% res = feGetRep(fe,'fiber res sig with mean',coords);
% res = feGetRep(fe,'fiber res sig with mean',voxelIndices);
%---------
% Residual signal fiber model prediction - demeaned measured signal
% with added mean signal. Res is returned as an array (nBvecs x nVoxel).
% This is used to reconstruct an image (volume) to be used for the
% refinement process.
% res = feGetRep(fe,'fiber res sig with mean voxel');
% res = feGetRep(fe,'fiber res sig with mean voxel',coords);
% res = feGetRep(fe,'fiber res sig with mean voxel',voxelIndices);
%---------
% Return a column vector of the proportion of variance explained in
% each voxel.
% R2byVox = feGetRep(fe,'voxr2');
% R2byVox = feGetRep(fe,'voxr2',coords);
%---------
% Return a column vector of the proportion of variance explained in
% each voxel. (Normalized to the squared mean diffusion signal in each voxel)
% R2byVox = feGetRep(fe,'voxr2zero');
% R2byVox = feGetRep(fe,'voxr2zero',coords);
%---------
% Return the percent of varince explained in each voxel.
% R2byVox = feGetRep(fe,'var exp by voxel');
% R2byVox = feGetRep(fe,'var exp by voxel',coords);
%---------
% Demeaned diffusion signal in each voxel.
% dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel');
% dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel',coords);
% dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel',vxIndex);
%---------
% Full (measured) signal in each voxel.
% dSigByVoxel = feGetRep(fe,'dSig full by Voxel');
% dSigByVoxel = feGetRep(fe,'dSig full by Voxel',coords);
% dSigByVoxel = feGetRep(fe,'dSig full by Voxel',vxIndex);
%---------
% Predicted signal by the full model in a set of voxeles.
% pSigByVoxel = feGetRep(fe, 'pSig full by voxel');
% pSigByVoxel = feGetRep(fe, 'pSig full by voxel',coords);
% pSigByVoxel = feGetRep(fe, 'pSig full by voxel',voxelIndex);
%---------
% A volume of RMSE values.
% RMSE = feGetRep(fe,'vox rmse')
% RMSE = feGetRep(fe,'vox rmse',coords)
% RMSE = feGetRep(fe,'vox rmse',vxIndex)
%---------
% A volume of RMSE values with a subset of fibers' weights set to 0.
% RMSE = feGetRep(fe,'vox rmse',fiberIndices)
% RMSE = feGetRep(fe,'vox rmse',fiberIndices,coords)
% RMSE = feGetRep(fe,'vox rmse',fiberIndices,voxelIndex)
%---------
% Fibers' residual signal by voxel.
% res = feGetRep(fe,'res sig fiber vox')
% res = feGetRep(fe,'res sig fiber vox',coords)
% res = feGetRep(fe,'res sig fiber vox',vxIndex)
%---------
% Full (measured) residual signal by voxel.
% res = feGetRep(fe,'res sig full vox')
% res = feGetRep(fe,'res sig full vox',coords)
% res = feGetRep(fe,'res sig full vox',vxIndex)
%---------
% Fibers' residual signal by voxel with mean signal added (with added
% isotropic component).
% This is used to compute the residual signal for the refinement.
% res = feGetRep(fe,'fiber res sig with mean vox')
% res = feGetRep(fe,'fiber res sig with mean vox',coords)
% res = feGetRep(fe,'fiber res sig with mean vox',vxIndex)
%---------
% Residual signal full model prediction from a multi-voxel fit.
% res = feGetRep(fe,'res sig full voxfit');
% res = feGetRep(fe, 'res sig full voxfit',coords);
% res = feGetRep(fe, 'res sig full voxfit',voxelIndex);
%---------
% Dimensions of the DW volume.
% dim = feGetRep(fe,'dims')
%---------
% Dimensions of the maps of parameters and results.
% dims = feGetRep(fe, 'mapsize')
%
% End of feGetRep.m parameters, 
% 
%
% Copyright (2013-2014), Franco Pestilli, Stanford University, pestillifranco@gmail.com.

val = [];

% Format the input parameters.
param = lower(strrep(param,' ',''));

% Start sorting the input and computing the output.
switch param
  case 'name'
    % Name of the current fe structure.
    %
    % name = feGetRep(fe,'name')
    val = fe.rep.name;
    
  case 'type'
    % The type of objes (always, fascicle evaluation)
    %
    % type = feGetRep(fe,'type')
    val = fe.rep.type; % Always fascicle evaluation
      
  case 'dwi'
    %  Load the diffusion weighted data.
    %
    % dwi = feGetRep(fe,'dwirepeat')
    val = dwiLoad(feGetRep(fe,'dwifile'));
        
  case 'dwifile'
    %  Load the repeated measure of the diffusion weighted data.
    %
    % dwiFile = feGetRep(fe,'dwirepeatfile')
    val = fe.path.dwifilerep;
    
  case {'bvecs'}
    % Diffusion directions.
    %
    % val = feGetRep(fe,'bvecs');
    val = fe.rep.bvecs;
    
  case {'bvecsindices'}
    % Indices to the diffusion directions in the DWi 4th Dimension.
    %
    % val = feGetRep(fe,'bvecs indices rep');
    val = fe.rep.bvecsindices;
    
  case {'nbvecs','nbvals'}
    % Number of B0's
    %
    % val = feGetRep(fe,'n bvals');
    val = length(feGetRep(fe,'bvals'));
    
  case {'bvals'}
    % B0 Values.
    %
    % bval = feGetRep(fe,'bvals')
    val = fe.rep.bvals;
    
  case {'diffusionsignalinvoxel','dsiinvox','dsigvox','dsigmeasuredvoxel'}
    % Returns a nVoxels X nBvecs array of measured diffusion signal
    %
    % val = feGetRep(fe,'dsiinvox');
    % val = feGetRep(fe,'dsiinvox',voxelsIndices);
    % val = feGetRep(fe,'dsiinvox',coords);
    val = fe.rep.diffusion_signal_img(feGet(fe,'voxelsindices',varargin),:)';
    
  case {'diffusionsignalinvoxeldemeaned','dsiinvoxdemeaned'}
    % Returns a nVoxels X nBvecs array of demeaned diffusion signal
    %
    % val =feGetRep(fe,'dsiinvoxdemeaned');
    % val =feGetRep(fe,'dsiinvoxdemeaned',voxelsIndices);
    % val =feGetRep(fe,'dsiinvoxdemeaned',coords);
    nBvecs       = feGetRep(fe,'nBvecs');
    voxelIndices = feGet(fe,'voxelsindices',varargin);
    val = fe.rep.diffusion_signal_img(voxelIndices,:) - repmat(mean(fe.rep.diffusion_signal_img(voxelIndices,:), 2),1,nBvecs);
    keyboard % THis seems to be worng
    
  case {'b0signalimage','b0vox'}
    % Get the diffusion signal at 0 diffusion weighting (B0) for this voxel
    %
    % val = feGetRep(fe,'b0signalimage');
    % val = feGetRep(fe,'b0signalimage',voxelIndex);
    % val = feGetRep(fe,'b0signalimage',coords);
    val = fe.rep.diffusion_S0_img(feGet(fe,'voxelsindices',varargin), :);
 
  case {'isoweights','weightsiso','meanvoxelsignal'}
    % Weights of the isotropic voxel signals, this is the mean signal in
    % each voxel.
    %
    % val = feGetRep(fe,'iso weights');
    % val = feGetRep(fe,'iso weights'coords)
    % val = feGetRep(fe,'iso weights',voxelIndices)
    val = feGet(fe,'Miso') \ feGetRep(fe,'dsig full')';
    val = val(feGet(fe,'voxelsindices',varargin));
     
  case {'dsigiso','isodsig'}
    % mean signalin each voxel, returned for each diffusion direction:
    % size(nBvecsxnVoxels, 1).
    %
    % val = feGetRep(fe,'iso disg');
    % val = feGetRep(fe,'iso disg',coords)
    % val = feGetRep(fe,'iso dsig',voxelIndices)
    val = feGetRep(fe,'iso weights');
    val = repmat(val,1,feGetRep(fe,'nbvecs'))';
    val = val(:);
    
  case {'dsigmeasured','dsigfull'}
    % Measured signal in VOI, this is the raw signal. not demeaned
    %
    % dSig = feGetRep(fe,'dSig full')
    val = fe.rep.diffusion_signal_img';
    val = val(:)';

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
    % dSig = feGetRep(fe,'dsigdemeaned');
    % dSig = feGetRep(fe,'dsigdemeaned',[1 10 100]);
    % dSig = feGetRep(fe,'dsigdemeaned',coords);
    nVoxels = feGet(fe,'nVoxels');
    nBvecs  = feGetRep(fe,'nBvecs');
    val     = (feGetRep(fe,'dsig measured') - reshape(repmat( ...
      mean(reshape(feGetRep(fe,'dsig measured'), nBvecs, nVoxels),1),...
      nBvecs,1), size(feGetRep(fe,'dsig measured'))))';
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
    % dSig = feGetRep(fe,'dsigrowssubset',voxelsList);
    dSig = feGetRep(fe,'dsigdemeaned');
    for vv = 1:length(voxelsList)
      val(feGet(fe,'voxel rows',vv)) = dSig(feGet(fe,'voxel rows',voxelList(vv)));
    end
    
  case {'totalr2'}
    % Return the global R2 (fraction of variance explained) of the full life
    % model.
    %
    % R2 = feGetRep(fe,'total r2');

    %     measured  = feGetRep(fe,'dsigdemeaned');
    %     predicted = feGetRep(fe,'fiber p sig');
    %     val = (1 - (sum((measured - predicted).^2 ) ./ ...
    %                 sum((measured - mean(measured)).^2) ));
    val = (1 - (sum((feGetRep(fe,'diffusion signal demeaned') - ...
      feGet(fe,'pSig fiber')).^2 ) ./ ...
      sum((feGetRep(fe,'diffusion signal demeaned') - ...
      mean(feGetRep(fe,'diffusion signal demeaned'))).^2) ));
  
  case {'totalr2voxelwise'}
    % Return the global R2 (fraction of variance explained) of the full life
    % model from a voxel-wise fit
    %
    % R2 = feGetRep(fe,'total r2');

    %     measured  = feGetRep(fe,'dsigdemeaned');
    %     predicted = feGetRep(fe,'p sig f voxel wise');
    %     val = (1 - (sum((measured - predicted).^2 ) ./ ...
    %                 sum((measured - mean(measured)).^2) ));
    val = (1 - (sum((feGetRep(fe,'diffusion signal demeaned')' - ...
      feGet(fe,'pSig f voxel wise')).^2 ) ./ ...
      sum((feGetRep(fe,'diffusion signal demeaned') - ...
      mean(feGetRep(fe,'diffusion signal demeaned'))).^2) ));
  
  case {'totpve'}
    % Total percent variance explained by data1 on data2
    %
    % R2 = feGetRep(fe,'tot pve data')
    val = 100 * feGetRep(fe,'total r2');

  case {'totalr2data'}
    % Return the global R2 (fraction of variance explained) of the full life
    % model in data set 2 by data set 1.
    %
    % R2 = feGetRep(fe,'total r2 data');

    %     measured  = feGetRep(fe,'dsigdemeaned');
    %     predicted = feGetRep(fe,'fiber p sig');
    %     val = (1 - (sum((measured - predicted).^2 ) ./ ...
    %                 sum((measured - mean(measured)).^2) ));
    val = (1 - (sum((feGetRep(fe,'diffusion signal demeaned') - ...
      feGet(fe,'diffusion signal demeaned')).^2 ) ./ ...
      sum((feGetRep(fe,'diffusion signal demeaned') - ...
      mean(feGetRep(fe,'diffusion signal demeaned'))).^2) ));
 
  case {'totpvedata'}
    % Total percent variance explained by data1 on data2
    %
    % R2 = feGetRep(fe,'tot pve data')
    val = 100 * feGetRep(fe,'total r2 data');
    
  case {'totalrmsedata'}
    % Global root mean squared error of data1 on data 2
    %
    % rmse = feGetRep(fe,'rmse data')
    val = sqrt(mean((feGetRep(fe,'diffusion signal demeaned') - ...
      feGet(fe,'diffusion signal demeaned')).^2));
      
  case {'totalrmseratio'}
    % Global ratio of root mean squared error of Model/Data
    %
    % rmse = feGetRep(fe,'rmse data ratio')
    val = feGetRep(fe,'total rmse') ./ feGetRep(fe,'total rmse data');
     
  case {'totalrmseratiovoxelwise'}
    % Global ratio of root mean squared error of Model/Data from a
    % voxel-wise fit
    %
    % rmse = feGetRep(fe,'totalrmseratiovoxelwise')
    val = feGetRep(fe,'total rmse voxelwise') ./ feGetRep(fe,'total rmse data');

  case {'totalrmse','rmsetotal'}
    % Root mean squared error of the LiFE fit to the whole data
    %
    % rmse = feGetRep(fe,'rmsetotal')
    val = sqrt(mean((feGetRep(fe,'diffusion signal demeaned') - ...
                     feGet(fe,'psigfiber')).^2));
  
  case {'totalrmsevoxelwise'}
    % Root mean squared error of the LiFE fit to the whole data from a
    % voxel-wise fit
    %
    % rmse = feGetRep(fe,'totalrmsevoxelwise')
    val = sqrt(mean((feGetRep(fe,'diffusion signal demeaned')' - ...
      feGet(fe,'psigfvoxelwise')).^2));

  case {'ressigfiber'}
    % Residual signal fiber prediction - measured_demeaned.
    %
    % res = feGetRep(fe,'res sig fiber')
    val = (feGetRep(fe,'dsigdemeaned')    - feGet(fe,'psig fiber'));
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val              = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'ressigfull'}
    % Residual signal full model prediction - measured signal.
    %
    % res = feGetRep(fe,'res sig full');
    % res = feGetRep(fe, 'res sig full',coords);
    % res = feGetRep(fe, 'res sig full',voxelIndex);
    val = (feGetRep(fe,'dsig full')    - feGet(fe,'psig full')');
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'ressigfullvoxelwise'}
    % Residual signal full model prediction from a multi-voxel fit.
    %
    % res = feGetRep(fe,'res sig full voxfit');
    % res = feGetRep(fe, 'res sig full voxfit',coords);
    % res = feGetRep(fe, 'res sig full voxfit',voxelIndex);
    %tic,val      = feGetRep(fe,'dsig demeaned') - feGet(fe,'psigfvoxelwise') + feGetRep(fe,'dsigiso');toc
    val      = feGetRep(fe,'dsig full') - feGet(fe,'psigfvoxelwise');
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
   case {'voxressigfullvoxelwise'}
    % Residual signal full model prediction from a multi-voxel fit.
    % ORganized in nvecsxn
    %
    % res = feGetRep(fe,'voxressigfullvoxelwise');
    % res = feGetRep(fe, 'voxressigfullvoxelwise',coords);
    % res = feGetRep(fe, 'voxressigfullvoxelwise',voxelIndex);
    %tic,val      = feGetRep(fe,'dsig demeaned') - feGet(fe,'psigfvoxelwise') + feGetRep(fe,'dsigiso');toc
    val = reshape(feGetRep(fe,'ressigfullvoxelwise'),feGetRep(fe,'nbvecs'),feGet(fe,'nvoxels'));

    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(:,feGet(fe,'voxelsindices',varargin));
    end
    
  case {'voxisoweights','voxweightsiso','voxmeanvoxelsignal'}
    % Weights of the isotropic voxel signals, this is the mean signal in
    % each voxel.
    % Organized nBvecs x nVoxels
    %
    % val = feGet(fe,'voxiso weights');
    % val = feGet(fe,'voxiso weights'coords)
    % val = feGet(fe,'voxiso weights',voxelIndices)
    val = feGet(fe,'Miso') \ feGetRep(fe,'dsig full')';
    val = repmat(val,1,feGet(fe,'nbvecs'))';

    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(:,feGet(fe,'voxelsindices',varargin));
    end
    
  case {'fiberressigwithmean'}
    % Residual signal fiber model prediction - demeaned measured signal
    % with added mean signal. VECTOR FORM.
    %
    % This is used to reconstruct an image (volume) to be used for the
    % refinement process.
    %
    % res = feGetRep(fe,'fiber res sig with mean');
    
    % predicted = feGet(fe,'psig fiber');
    % measured  = feGetRep(fe,'dsigdemeaned');
    % val = (measured_full - predicted demeaned);
    val = (feGetRep(fe,'dsig full')' - feGet(fe,'psig fiber')); 
    
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end
    
  case {'predictedfullsignalvoxel','psigfullvox'}
    % Predicted signal by the full model in a set of voxeles.
    %
    % Vox subfield stores per voxel within the VOI
    % The pSig has size of the dwi.
    % It is stored as val = pSig(X,Y,Z,Theta)
    %
    % pSig = feGetRep(fe, 'pSig full by voxel');
    % pSig = feGetRep(fe, 'pSig full by voxel',coords);
    % pSig = feGetRep(fe, 'pSig full by voxel',voxelIndex);
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'pSig full'), nBvecs, nVoxels);
      
  case {'voxelr2','r2vox','voxr2','r2byvoxel'}
    % Return a column vector of the proportion of variance explained in
    % each voxel.
    %
    % R2byVox = feGetRep(fe,'voxr2');
    % R2byVox = feGetRep(fe,'voxr2',coords);
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    nBvecs    = feGetRep(fe,'nBvecs');
    val       = (1 - (sum((measured - predicted).^2 ) ./ ...
      sum((measured - repmat(mean(measured), nBvecs,1)).^2 ) ));
    val       = val(feGet(fe,'return voxel indices',varargin));
 
  case {'voxelss','ssvox','voxss','sumofsquaresbyvoxel'}
    % Return a column vector of sum of squares error of the model 
    % prediction in each voxel
    %
    % ssem = feGetRep(fe,'voxss');
    % ssem = feGetRep(fe,'voxss',coords);
    predicted = feGet(fe,'pSig f vox');    
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    val       = sum((measured - predicted).^2 );
    val       = val(feGet(fe,'return voxel indices',varargin));
 
  case {'voxelssdata','ssvoxdata','voxssdata','sumofsquaresbyvoxeldata'}
    % Return a column vector of sum of squares error of the data in
    % each voxel
    %
    % ssed = feGetRep(fe,'voxssedata');
    % ssed = feGetRep(fe,'voxssdata',coords);
    measured1 = feGet(fe,'dSig demeaned by voxel');
    measured2  = feGetRep(fe,'dSig demeaned by voxel');
    val       = sum((measured2 - measured1).^2 );
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelr2data','r2voxdata','voxr2data','r2byvoxeldata'}
    % Return a column vector of the proportion of variance explained in
    % each voxel in data set 2 given the data in dataset 1.
    %
    % R2byVox = feGetRep(fe,'voxr2data');
    % R2byVox = feGetRep(fe,'voxr2data',coords);
    measured1  = feGet(fe,'dSig demeaned by voxel');
    measured2 = feGetRep(fe,'dSig demeaned by voxel');
    nBvecs    = feGetRep(fe,'nBvecs');
    val       = (1 - (sum((measured2 - measured1).^2 ) ./ ...
      sum((measured2 - repmat(mean(measured2), nBvecs,1)).^2 ) ));
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelr2zero','r2voxzero','voxr2zero','r2byvoxelzero'}
    % Return a column vector of the proportion of variance explained in
    % each voxel. (Normalized to the squared diffusion signal in each voxel)
    %
    % R2byVox = feGetRep(fe,'voxr2zero');
    % R2byVox = feGetRep(fe,'voxr2zero',coords);
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    val       = (1 - (sum((measured - predicted).^2 ) ./ sum(measured.^2)));
    val       = val(feGet(fe,'return voxel indices',varargin));
     
  case {'voxelr2pearson','r2voxpearson','voxr2corr','r2byvoxelpearson'}
    % Return a column vector of the squre of the pearson correlation coefficient
    %
    % R2byVox = feGetRep(fe,'voxr2corr');
    % R2byVox = feGetRep(fe,'voxr2corr',coords);
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    val       = corrcoef(measured - predicted).^2;
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelr2zerodata','r2voxzerodata','voxr2zerodata','r2byvoxelzerodata'}
    % Return a column vector of the proportion of variance explained in
    % each voxel in dataset 1 given the data in data set 2. 
    % (Normalized to the squared diffusion signal in each voxel)
    %
    % R2byVox = feGetRep(fe,'voxr2zerodata');
    % R2byVox = feGetRep(fe,'voxr2zerodata',coords);
    measured  = feGet(fe,'dSig demeaned by voxel');
    measured2 = feGetRep(fe,'dSig demeaned by voxel');
    val       = (1 - (sum((measured - measured2).^2 ) ./ sum(measured.^2)));
    val       = val(feGet(fe,'return voxel indices',varargin));

  
  case {'voxelr2voxelwise','r2voxvoxelwise','voxr2voxelwise','r2voxelwisebyvoxel'}
    % Return a column vector of the proportion of variance explained in
    % each voxel.
    %
    % R2byVox = feGetRep(fe,'voxr2voxelwise');
    % R2byVox = feGetRep(fe,'voxr2voxelwise',coords);
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f voxel wise by voxel');
    nBvecs    = feGetRep(fe,'nBvecs');
    val       = (1 - (sum((measured - predicted).^2 ) ./ ...
      sum((measured - repmat(mean(measured), nBvecs,1)).^2 ) ));
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelr2zerovoxelwise','r2zerovoxelwisebyvoxel','voxr2zerovoxelwise'}
    % Return a column vector of the proportion of variance explained in
    % each voxel. (Normalized to the squared diffusion signal in each voxel)
    %
    % R2byVox = feGetRep(fe,'voxr2zerovoxelwise');
    % R2byVox = feGetRep(fe,'voxr2zerovoxelwise',coords);
    measured  = feGetRep(fe,'dSig demeaned by voxel');
    predicted = feGet(fe,'pSig f voxel wise by voxel');
    val       = (1 - (sum((measured - predicted).^2 ) ./ sum(measured.^2)));
    val       = val(feGet(fe,'return voxel indices',varargin));    
    
  case {'voxelvarianceexplained','varexpvox','voxvarexp','varexpbyvoxel'}
    % Return the percent of varince explained in each voxel.
    %
    % R2byVox = feGetRep(fe,'var exp by voxel');
    % R2byVox = feGetRep(fe,'var exp by voxel',coords);
    val = 100.*feGetRep(fe,'voxr2');
    val = val(feGet(fe,'return voxel indices',varargin));
   
  case {'voxelvarianceexplaineddata','varexpvoxdata','voxvarexpdata','varexpbyvoxeldata'}
    % Return the percent of variance explained in each voxel in data set 1
    % given the data in data set 2.
    %
    % R2byVox = feGetRep(fe,'var exp by voxel data');
    % R2byVox = feGetRep(fe,'var exp by voxel data',coords);
    val = 100.*feGetRep(fe,'voxr2 data');
    val = val(feGet(fe,'return voxel indices',varargin));
 
  case {'dsigdemeanedbyvoxel','dsigdemeanedvox'}
    % Demeaned diffusion signal in each voxel
    %
    % dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel');
    % dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel',coords);
    % dSigByVoxel = feGetRep(fe,'dsigdemeaned by Voxel',vxIndex);
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'dsigdemeaned'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin));
    
  case {'dsigfullbyvoxel','dsigfullvox'}
    % Full (measured) signal in each voxel
    %
    % dSigByVoxel = feGetRep(fe,'dSig full by Voxel');
    % dSigByVoxel = feGetRep(fe,'dSig full by Voxel',coords);
    % dSigByVoxel = feGetRep(fe,'dSig full by Voxel',vxIndex);
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'dSig full'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
    
  case {'voxelrmse','voxrmse'}
    % A volume of RMSE values
    %
    % RMSE = feGetRep(fe,'vox rmse')
    % RMSE = feGetRep(fe,'vox rmse',coords)
    % RMSE = feGetRep(fe,'vox rmse',vxIndex)
    measured  = feGetRep(fe,'dsigdemeaned by voxel');
    predicted = feGet(fe,'pSig f vox');
    val       = sqrt(mean((measured - predicted).^2,1));
    val       = val(feGet(fe,'return voxel indices',varargin));
    
  case {'voxelrmsevoxelwise','voxrmsevoxelwise'}
    % A volume of RMSE values optained with the voxel-wise (vw) fit.
    %
    % RMSE = feGetRep(fe,'vox rmse vw')
    % RMSE = feGetRep(fe,'vox rmse vw',coords)
    % RMSE = feGetRep(fe,'vox rmse vw',vxIndex)
    measured  = feGetRep(fe,'dsigdemeaned by voxel');
    predicted = feGet(fe,'pSig f voxel wise by voxel');
    val       = sqrt(mean((measured - predicted).^2,1));
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelrmseratio','voxrmseratio'}
    % A volume with the ratio of the RMSE of model/data
    %
    % rmseRatio = feGetRep(fe,'vox rmse ratio')
    % rmseRatio = feGetRep(fe,'vox rmse ratio',coords)
    % rmseRatio  = feGetRep(fe,'vox rmse ratio',vxIndex)
    rmseData  = feGetRep(fe,'vox rmse data');
    rmseModel = feGetRep(fe,'vox rmse');
    val       = rmseModel ./ rmseData;
    val       = val(feGet(fe,'return voxel indices',varargin));
  
  case {'prmseratio','proportionrmseratio'}      
    % The probability of a ratio-value in the volume.
    % Default across 25 log-distributed bins between [.5,2]
    %
    % rmseRatio = feGetRep(fe,'p rmse ratio')    
    %
    % % Change the bins over whih the proportions are computed:
    % bins      = logspace(log10(.25),log10(4),50)
    % rmseRatio = feGetRep(fe,'p rmse ratio',bins)
    if isempty(varargin)
        bins = logspace(log10(.5),log10(2),25);
    end
    % Extract the rmse ratio in each voxel
    Rrmse  = feGetRep(fe,'vox rmse ratio');
    % Compute the number of occurrences for a range of values.
    [val(1,:),val(2,:)]  = hist(Rrmse,bins);
    val(1,:)             = val(1,:)./sum(val(1,:));
    
  case {'voxelrmseratiovoxelwise','voxrmseratiovoxelwise'}
    % A volume with the ratio of the RMSE of model/data
    %
    % rmseRatio = feGetRep(fe,'vox rmse ratio voxel wise')
    % rmseRatio = feGetRep(fe,'vox rmse ratio voxel wise',coords)
    % rmseRatio  = feGetRep(fe,'vox rmse ratio voxel wise',vxIndex)
    rmseData  = feGetRep(fe,'vox rmse data');
    rmseModel = feGetRep(fe,'vox rmse voxel wise');
    val       = rmseModel ./ rmseData;
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelrmsedata','voxrmsedata'}
    % A volume of RMSE values from data set 1 to data set 2
    %
    % RMSE = feGetRep(fe,'vox rmse data')
    % RMSE = feGetRep(fe,'vox rmse data',coords)
    % RMSE = feGetRep(fe,'vox rmse data',vxIndex)
    measured  =    feGet(fe,'dsigdemeaned by voxel');
    measured2 = feGetRep(fe,'dsigdemeaned by voxel');
    val       = sqrt(mean((measured - measured2).^2,1));
    val       = val(feGet(fe,'return voxel indices',varargin));

  case {'voxelrmsetest','voxrmsetest'}
    % A volume of RMSE values with a subset of fibers' weights set to 0.
    %
    % RMSE = feGetRep(fe,'vox rmse',fiberIndices)
    % RMSE = feGetRep(fe,'vox rmse',fiberIndices,coords)
    % RMSE = feGetRep(fe,'vox rmse',fiberIndices,voxelIndex)
    measured  = feGetRep(fe,'dsigdemeaned by voxel');
    % Reshape the predicted signal by voxles
    predicted = reshape(feGet(fe,'pSig fiber test voxel wise',varargin{1}), feGetRep(fe,'nBvecs'), feGet(fe,'nVoxels'));   
    val       = sqrt(mean((measured - predicted).^2,1));
    
    if length(varargin) == 2
      val       = val(feGet(fe,'return voxel indices',varargin));
    end
    
  case {'voxelrmsetestvoxelwise','voxrmsetestvoxelwise'}
    % A volume of RMSE values with a subset of fibers' weights set to 0.
    %
    % RMSE = feGetRep(fe,'voxel rmse test voxel wise',fiberIndices)
    % RMSE = feGetRep(fe,'voxel rmse test voxel wise',fiberIndices,coords)
    % RMSE = feGetRep(fe,'voxel rmse test voxel wise',fiberIndices,voxelIndex)
    measured  = feGetRep(fe,'dsigdemeaned by voxel');
    % Reshape the predicted signal by voxles
    predicted = reshape(feGet(fe,'pSig fiber test voxel wise',varargin{1}), feGetRep(fe,'nBvecs'), feGet(fe,'nVoxels'));   
    val       = sqrt(mean((measured - predicted).^2,1));
    
    if length(varargin) == 2
      val       = val(feGet(fe,'return voxel indices',varargin));
    end
    
  case {'residualsignalfibervoxel','resfibervox'}
    % Fibers' residual signal by voxel
    %
    % res = feGetRep(fe,'res sig fiber vox')
    % res = feGetRep(fe,'res sig fiber vox',coords)
    % res = feGetRep(fe,'res sig fiber vox',vxIndex)
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'res sig fiber'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
    
  case {'residualsignalfullvoxel','resfullvox'}
    % Full (measured) residual signal by voxel
    %
    % res = feGetRep(fe,'res sig full vox')
    % res = feGetRep(fe,'res sig full vox',coords)
    % res = feGetRep(fe,'res sig full vox',vxIndex)
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'res sig full'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))';
      
  case {'fiberressigwithmeanvoxel'}
    % Residual signal fiber model prediction - demeaned measured signal
    % with added mean signal. VOXEL FORM (nBvecs x nVoxel).
    %
    % This is used to reconstruct an image (volume) to be used for the
    % refinement process.
    %
    % res = feGetRep(fe,'fiber res sig with mean');
    
    % predicted = feGet(fe,'psig fiber');
    % measured  = feGetRep(fe,'dsigdemeaned');
    % val = (measured     - predicted) + (measured_full - measured);
    val     = feGetRep(fe,'fiberressigwithmean');
    nBvecs  = feGetRep(fe,'nBvecs');
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
    % res = feGetRep(fe,'fiber res sig with mean vox')
    % res = feGetRep(fe,'fiber res sig with mean vox',coords)
    % res = feGetRep(fe,'fiber res sig with mean vox',vxIndex)
    nBvecs  = feGetRep(fe,'nBvecs');
    nVoxels = feGet(fe,'n voxels');
    val     = reshape(feGetRep(fe,'fiber res sig with mean'), nBvecs, nVoxels);
    val     = val(:,feGet(fe,'return voxel indices',varargin))'; 

  case {'volumesize','dims','dim','imagedim'}
    % Dimensions of the DW volume.
    %
    % dim = feGetRep(fe,'dims')
    val = fe.rep.imagedim;
  case {'mapsize'}
    % Dimensions of the maps of parameters and results.
    %
    % dims = feGetRep(fe, 'mapsize')
    val = fe.rep.imagedim(1:3);
    
  otherwise
    help('feGetRep')
    fprintf('[feGetRep] Unknown parameter << %s >>...\n',param);
    keyboard
end

end % END MAIN FUNCTION
