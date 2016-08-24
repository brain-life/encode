function fe = feCreate
% Create a linear fascicle evaluation structure
%
%  fe = feCreate
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Book-keeping
fe.name  = 'default'; % This one's name
fe.type  = 'faseval'; % Always fascicle evaluation

% Overview
fe.life    = [];  % Structure of parameters and results from LIFE analysis
fe.fg      = [];  % Fiber group candidate fascicles, uses fgSet/Get
fe.roi     = [];  % Cell array of regions of interest where we evaluate

% Path to files
fe.path.dwifile = [];  % Diffusion weighted file used for testing results
fe.path.dtfile  = [];  % Diffusion weighted file used for testing results
fe.path.savedir = [];  % Top directory under which all files will be saved
fe.path.anatomy = [];  % 3D high-resolution anatomical file
fe.path.dwifilerep = [];  % File used for cross-validating the results. 
                          % A secodn measuremnt ideally acquried in the 
                          % same session witht the same scanning
                          % parameters.

% Life Algorithm parameters
%fe.life.Mfiber       = []; % Fiber portion of A matrix
%fe.life.voxel2FNpair = []; % Voxels to fiber/node pair, fefgGet(fg,'voxel 2 fiber node pair')
%fe.life.M.orient          = []; % Dictionary of orientations
%fe.life.M.DictSig         = []; % Dictionary of Demeaned Signals
fe.life.M.Phi             = []; % Indication Tensor: sparse-(nFibers x nAtoms x nVoxels) 3way array
%fe.life.dSig         = []; % Signal of fibers alone (demeaned).
%fe.life.Mfiber       = [];
%fe.life.voxel2FNpair = {};
%fe.life.dSig         = [];
fe.life.xform        = [];
fe.life.fibers       = [];
fe.life.diffusion_signal_img = [];
fe.life.diffusion_S0_img     = [];
fe.life.bvecs                = [];
fe.life.bvals                = [];
fe.life.bvecsindices         = [];
fe.life.imagedim             = [];
fe.life.xform          = [];  % Transforms between coords for fg, roi, and dwi data
fe.life.xform.img2acpc = [];
fe.life.xform.acpc2img = [];

% Repetition file for cross-validation
fe.rep = [];

return
