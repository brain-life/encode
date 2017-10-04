function fe = feSet(fe,param,val,varargin)
% Set fascicle evaluation parameters.
%
%   fe = feSet(fe,param,val,varargin)
%
%----------
% feSet(fe,'bvecsindices');
%----------   
% Set the a second diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'diffusion signal repeat',image_vals);
%----------
% Set the the s0 (no diffusion direction measurement) for a second
% diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'b0signalrepeat',image_vals);
%----------
%
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Check for input parameters
if notDefined('fe'),    error('fe structure required'); end
if notDefined('param'), error('param required'); end
if ~exist('val','var'), error('Value required'); end

% Squeeze out spaces and force lower case
param = mrvParamFormat(param);

%%
switch param
  % Book-keeping
  case 'name'
    fe.name  = val;        % This one's name
  case 'type'
    fe.type  = 'faseval'; % Always fascicle evaluation
  case 'savedir'
    fe.path.savedir  = val; % Always fascicle evaluation
  
    % Set top level structure, not just single slot
  case 'life'
    fe.life  = val;  % Structure of parameters and results from LIFE analysis
    
  case 'fgfromacpc'
    % Fiber group candidate fascicles, Connectome.
    % Everything isin img coordinates in LiFE
    fe.fg  = dtiXformFiberCoords(val, fe.life.xform.acpc2img,'img');
    % Must clear when we change the fg
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);
    
  case {'fgimg', 'fg'}
    % Fiber group candidate fascicles, Connectome.
    % Everything is in img coordinates in LiFE
    %
    fe.fg  = val;
    % Must clear when we change the fg
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);

  case {'fgtensors','tensors','fgq','q'}
    % fe = feSet(fe, 'tensors', tensors);    - set the passed tensors
    fe.life.fibers.tensors = val;

  case 'roi'
    % Cell array of regions of interest where we evaluate
    % Must clear v2fnp
    fe.roi   = val;
    fe       = feSet(fe,'v2fnp',[]);
    
  case {'roifromfg','fgroi','roifg'}
    name   = sprintf('roi_%s',fe.fg.name); 
    randColor = rand(1,3);
    fe.roi = dtiNewRoi(name,randColor,fefgGet(feGet(fe,'fg img'),'unique image coords'));
    
  case 'xform'
    fe.life.xform = val;  % Transforms between coords for fg, roi, and dwi data
     
  %% Diffusion data related parameters
  case {'bvecs','diffusionbvecs'}
    % feSet(fe,'bvecs');
    fe.life.bvecs = val;
  case {'bvecsindices','diffusionimagesindicesindwivolume'}
    % feSet(fe,'bvecsindices');
    fe.life.bvecsindices = val;
  case {'bvals','diffusionbvals'}
    fe.life.bvals = val;    
  case {'diffusionsignalimage','dsi', 'diffusion_signal_img'}
    fe.life.diffusion_signal_img = val;
  case {'b0signalimage','b0img', 'diffusion_s0_im','s0image'}
    fe.life.diffusion_S0_img = val;
  case {'usedvoxels'}
    fe.life.usedVoxels = val;
  case {'modeltensor'}
    fe.life.modelTensor = val;
  case {'roivoxels','roicoords'}
    % What space?  What form for the coords?
    % Always in IMG coords in LiFE.
    fe.roi.coords = val;
  
    
    %% The LiFE model
  case 'mfiber'
    fe.life.Mfiber = val;             % Fiber portion of M matrix
  case {'measuredsignalfull', 'dsigmeasured'}      % Measured signal in ROI
    fe.life.dSig  = val;
  case 'fit'
    fe.life.fit = val;
  case 'voxfit'
    fe.life.voxfit = val;
  case 'xvalfit'
    fe.life.xvalfit = val;
  
    %% Connectome fibers information.
  case {'numberofuniquefibersineachvoxel','uniquefibersnum','numberofuniquefibers','numuniquef'}
    fe.life.fibers.unique.num = val;
  case {'indextouniquefibersineachvoxel','uniquefibersindex','uniqueindex','indexesofuniquefibers','indexuniquef','uniquefibers'}
    fe.life.fibers.unique.index = val;
  case {'numberoftotalfibersineachvoxel','totalfibernmber','fibersnum','numberoffibers','numf','numfibers'}
    fe.life.fibers.total.num = val;
  case {'indexoftotalfibersineachvoxel','totalfiberindex','fibersbyvox','fibersinvox'}
    fe.life.fibers.total.index = val;
  case {'voxel2fibernodepairs','v2fnp'}
    % This has to be cleared whenever we change fg or roi
    fe.life.voxel2FNpair = val;
    % Spatial coordinate transforms for voxels and fg to coordinate frames
  case {'xformimg2acpc','img2acpc','img2acpcxform'}
    fe.life.xform.img2acpc = val;
  case {'xformacpc2img','acpc2img','acpc2imgxform'}
    fe.life.xform.acpc2img = val;
  case {'size','imgsize','volumesize','dims','dim'}
    fe.life.imagedim = val;
    
    %% Diffusion data reapeted measure parameters
  case 'dwirepeatfile'
    fe.path.dwifilerep = val;  % Diffusion weighted file used for testing results
  case {'diffusionsignalimagerepeat'}
    % Set the a second diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'diffusion signal repeat',image_vals);
    fe.rep.diffusion_signal_img = val;
  case {'s0imagerepeat'}
    % Set the the s0 (no diffusion direction measurement) for a second
    % diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'b0signalrepeat',image_vals);
    fe.rep.diffusion_S0_img = val;
  case {'bvecsrepeat','diffusionbvecsrepeat'}
    % feSet(fe,'bvecsrepeat');
    fe.rep.bvecs = val;
  case {'bvecsindicesrepeat','diffusionimagesindicesindwivolumerepeat'}
    % feSet(fe,'bvecsindicesrepeat');
    fe.rep.bvecsindices = val;
  case {'bvalsrepeat','diffusionbvalsrepeat'}
    % fe = feSet(fe,'bvalsrepeat')
    fe.rep.bvals = val;
  case {'imgsizerepeat'}
    fe.rep.imagedim = val;
    
  case {'anatomyfile'}
    fe.path.anatomy = val;
  case 'dwifile'
    fe.path.dwifile = val;  % Diffusion weighted file used for testing results
  case 'dtfile'
    fe.path.dtfile = val;  % Diffusion weighted file used for testing results
  case 'dictionaryparameters'
    fe.life.M.Nphi = val{1};
    fe.life.M.Ntheta = val{2};
    fe.life.M.orient = val{3};
    fe.life.M.Dict = val{4};   
    fe.life.M.DictSig = val{5};   
  case 'gradients'
    fe.life.fibers.grad = val{1};
    fe.life.fibers.Nelem = val{2};
  case 'indicationtensor'
    fe.life.M.Phi = val;
  case 'atomslist'   %% To revise if can be avoided
    fe.life.fibers.atoms_list = val;
  case 'voxelslist'  %% To revise if can be avoided
    fe.life.fibers.voxels_list = val;
  case 'nelem'   %% To revise if can be avoided
    fe.life.fibers.Nelem = val;
  case 'ms0'
    fe.life.M.S0 = val;
  case 'curvature'
    fe.fg.Cur = val{1};
    fe.fg.Indication = val{2};
  case 'torsion'
    fe.fg.Tor = val{1};  
    fe.fg.Indication = val{2}; 
  case 'tracts_info'
    Ntracts = size(val.names,2);
    for n=1:Ntracts
        fe.life.M.tracts{n}.ind = find(val.index==n);
        fe.life.M.tracts{n}.name = val.names{n};
    end
    fe.life.M.tracts{Ntracts+1}.ind = find(val.index==0);
    fe.life.M.tracts{Ntracts+1}.name = 'not a tract'; 
            
  otherwise
    error('Unknown parameter %s\n',param);
end

end
