function fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeated,anatomyFile,varargin)
% Initialize a new connectome (fe) structure. 
%
%    fe = feConnectomeInit(dwiFile,dtFile,fgFileName,feFileName,savedir,dwiFileRepeated,anatomyFile,varargin);
%    
% We allow a set of (paramName,val) pairs in the varargin that will be
% executed as fe = feSet(fe,paramName,val)
%
% Example: 
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

%feOpenLocalCluster

% Intialize the fe structure.
fe = feCreate;

% Set the based dir for fe, this dire will be used 

if notDefined('savedir')
    if isstruct(fgFileName)
        savedir = fullfile(fileparts(fgFileName.name),'life');
    else
        savedir = fullfile(fileparts(fgFileName),'life');
    end
end
fe = feSet(fe,'savedir',savedir);

% Set the xforms (transformations from diffusion data to acpc)
tempNi = niftiRead(dwiFile);
fe = feSet(fe, 'img2acpc xform', tempNi.qto_xyz);
fe = feSet(fe, 'acpc2img xform', inv(tempNi.qto_xyz));
clear tempNi

% Set up the fe name
if isstruct(fgFileName),  n  = fgFileName.name;
else                   [~,n] = fileparts(fgFileName);
end

if notDefined('feFileName'),
  feFileName = sprintf('%s-%s', datestr(now,30),n);
end
fe = feSet(fe, 'name',feFileName);

% Load a connectome
if isstruct(fgFileName), fg = fgFileName; clear fgFileName
else % A file name was passed load the fibers from disk
  fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgFileName)
  fg = fgRead(fgFileName);
end

% Set fg in the fe structure identifying the fg coordinate frame.
% Everything in LiFE is in img coordinates, but everyting in mrDiffusion is in acpc.  
% So here we assume the fibers are read in acpc and we xform them in img.
fe = feSet(fe,'fg from acpc',fg);

% Precompute the canonical tensors for each node in each fiber (NOT NEEDED ANYMORE).

% Set model for canonical tensor
%tic
if ~isempty(varargin)
  N = varargin{1}(1);
  axialDiffusion  = varargin{2}(1);
  radialDiffusion = varargin{2}(2);  
  Compute_matrix_M = varargin{3};
else % Default to stick and ball
  N = 360;   
  axialDiffusion  = 1;
  radialDiffusion = 0;
  Compute_matrix_M = 0;
end
dParms(1) =  axialDiffusion; 
dParms(2) = radialDiffusion; 
dParms(3) = radialDiffusion;
Nphi = N;
Ntheta = N;

fe = feSet(fe,'model tensor',dParms);

% NOT NEEDED ANYMORE
%fprintf('\n[%s] Computing fibers'' tensors... ',mfilename); 
%fe = feSet(fe, 'tensors', feComputeCanonicalDiffusion(fe.fg.fibers, dParms));  
%toc

%% This should be changed (computation of ROI)
fe = feSet(fe,'roi fg',[]);
clear fg

%% The following is not needed anymore
% When the ROI is set to empty, LiFE uses all of the voxels within the
% connectome as the ROI.

% % NEW: Compuet direction of fibers (gradient) and Maxelem (Maximum number of nodes per fiber)
% tic
% fprintf('\n[%s] Computing fiber directions (gradient) ...',mfilename); 
% fe = feSet(fe, 'gradients', feComputeGradients(fe.fg.fibers));  
% fprintf('took: %2.3fs.\n',toc)
% 
% % We disregard fibers that have identical trajectories within the ROI.
% roi = feGet(fe,'roi coords');
% fe  = feSet(fe,'voxel 2 fiber node pairs',fefgGet(feGet(fe,'fg img'),'v2fn',roi));
% clear roi;
% 
% % This compute the unique fibers per voxels
% fe  = feGetConnectomeInfo(fe);

% Install the information about the diffusion data.
fe = feConnectomeSetDwi(fe,dwiFile,0);

% Install the information about a repeated measurement of the diffusion
% data.
if ~notDefined('dwiFileRepeated')
  fe = feConnectomeSetDwi(fe,dwiFileRepeated,1);
end

% Install the path tot he anatomical high-resolution file.
if ~notDefined('anatomyFile')
  fe = feSet(fe,'anatomy file',anatomyFile);
end

%% NEW: Precompute Dictionary of orientations and canonical demeaned signals
% Define discretization steps for building Dictionary of canonical difussivities
% These numbers represents the steps in which the range [0,pi] is divided
% for spherical coordinates phi: azimuth and theta: pi/2-elevetion

fe = feBuildDictionaries(fe,Nphi,Ntheta);

%% NEW: The previous very large matrix M was replaced by a sparse multiway decomposition
% Build LiFE tensors and key connection matrices
fe = feConnectomeBuildModel(fe,Compute_matrix_M);


fprintf(['\n[%s] fe-structure Memory Storage:',ByteSize(fe),'\n'],mfilename);

return
