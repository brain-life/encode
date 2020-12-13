function fe = feConnectomeInit(dwiFile,fgFileName,dwiFileRepeated,anatomyFile,N,tensor,kurtosis,reset)
% Initialize a new connectome (fe) structure. 
%
%    fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeated,anatomyFile,varargin);
%    
% We allow a set of (paramName,val) pairs in the varargin that will be
% executed as fe = feSet(fe,paramName,val)
%
%  Copyright (2020) Indiana University
%
%  Franco Pestilli frakkopesto@gmail.com and 
%  Cesar F. Caiafa ccaiafa@gmail.com

%feOpenLocalCluster

if(~exist('N', 'var') || isempty(N))
    N = 360;
end

if(~exist('tensor', 'var') || isempty(tensor))
    axialDiffusion = 1; radialDiffusion = 0;
    tensor = [ axialDiffusion radialDiffusion radialDiffusion ];
end

if(~exist('kurtosis', 'var') || isempty(kurtosis))
    kurtosis = zeros(15,1);
end

if(~exist('reset', 'var') || isempty(reset))
    reset = true;
end

% Intialize the fe structure.
fe = feCreate;

% Set the based dir for fe, this dir will be used 
% if notDefined('savedir')
%     if isstruct(fgFileName)
%         savedir = fullfile(fileparts(fgFileName.name),'life');
%     else
%         savedir = fullfile(fileparts(fgFileName),'life');
%     end
% end
% fe = feSet(fe,'savedir',savedir);

% Set the xforms (transformations from diffusion data to acpc)
tempNi = niftiRead(dwiFile);
fe = feSet(fe, 'img2acpc xform', tempNi.qto_xyz);
fe = feSet(fe, 'acpc2img xform', inv(tempNi.qto_xyz));

% set the offset for files as well, then clear it
%fe = feSet(fe, 'dwi offset', tempNi.qto_xyz(1:3,4)');
clear tempNi

% % Set up the fe name
% if isstruct(fgFileName),  n  = fgFileName.name;
% else                   [~,n] = fileparts(fgFileName);
% end
%
% if notDefined('feFileName')
%   feFileName = sprintf('%s-%s', datestr(now,30),n);
% end
% fe = feSet(fe, 'name',feFileName);

%% Load a connectome from disk
if isstruct(fgFileName), fg = fgFileName; clear fgFileName
else % A file name was passed load the fibers from disk
  fprintf('\n[%s]\n loading fiber from file: %s\n',mfilename,fgFileName)
  fg = fgRead(fgFileName);
end
fe = feSet(fe,'fg from acpc',fg);
fe = feSet(fe,'roi fg',[]); clear fg

%% Diffusion data (Y)
% Install the information about the diffusion data.
fe = feConnectomeSetDwi(fe,dwiFile,0);

% Information about a repeated measurement of the diffusion data.
if ~notDefined('dwiFileRepeated')
  fe = feConnectomeSetDwi(fe,dwiFileRepeated,1);
end

%% Anatomy Install the path to the anatomical high-resolution file.
if ~notDefined('anatomyFile')
    
    fe = feSet(fe,'anatomy file',anatomyFile);
    
    try
        tempNi = niftiRead(anatomyFile);
        fe = feSet(fe, 'anat2acpcxform', tempNi.qto_xyz);
        fe = feSet(fe, 'acpc2anatxform', inv(tempNi.qto_xyz));
        clear tempNi
    catch
        error('Anatomy file unable to be opened with niftiRead()');
    end
       
    % add check to file headers for consistent orientation between dwi / anat
    
    % load minimum data for check - ignore repeat dwi for now
    doff = feGet(fe, 'dwi offset');
    aoff = feGet(fe, 'anat offset');
    
    % check x dimension is consistently signed
    if ~(all([ doff(1) aoff(1) ] < 0) || all([ doff(1) aoff(1) ] > 0))
        if reset
            warning('POTENTIAL X-FLIP - bvecs flipped on x-axis to correct.');
            bvecs = feGet(fe, 'bvecs');
            bvecs(:,1) = -bvecs(:,1);
            fe = feSet(fe, 'bvecs', bvecs);
        else
            warning('POTENTIAL X-FLIP - qoffset_x is different between the dwi and anatomy files.');
        end
    end
    
    % check y dimension is consistently signed
    if ~(all([ doff(2) aoff(2) ] < 0) || all([ doff(2) aoff(2) ] > 0))
        if reset
            warning('POTENTIAL Y-FLIP - bvecs flipped on y-axis to correct.');
            bvecs = feGet(fe, 'bvecs');
            bvecs(:,2) = -bvecs(:,2);
            fe = feSet(fe, 'bvecs', bvecs);
        else
            warning('POTENTIAL Y-FLIP - qoffset_y is different between the dwi and anatomy files.');
        end
    end
    
    % check z dimension is consistently signed
    if ~(all([ doff(3) aoff(3) ] < 0) || all([ doff(3) aoff(3) ] > 0))
        if reset
            warning('POTENTIAL Z-FLIP - bvecs flipped on z-axis to correct.');
            bvecs = feGet(fe, 'bvecs');
            bvecs(:,3) = -bvecs(:,3);
            fe = feSet(fe, 'bvecs', bvecs);
        else
            warning('POTENTIAL Z-FLIP - qoffset_z is different between the dwi and anatomy files.');
        end
    end
    
end

%% check / fix the tensor model for multishell input

nshells = feGet(fe, 'nshells');
if (size(tensor, 1) == 1) && nshells > 1
    tensorParams = repmat(tensor, [ nshells 1 ]);
else
    tensorParams = tensor;
end

% if data is single shell, zero out kurtosis model regardless of user input
if nshells == 1
    kurtosisParams = nan(15,1);
else
    kurtosisParams = kurtosis;
end

%% Precompute Dictionary of orientations and canonical demeaned signals
% Define discretization steps for building Dictionary of canonical difussivities
% These numbers represents the steps in which the range [0,pi] is divided
% for spherical coordinates phi: azimuth and theta: pi/2-elevetion

Nphi = N;
Ntheta = N;
fe = feSet(fe,'model tensor',tensorParams);
fe = feSet(fe,'model kurtosis',kurtosisParams);
% if single shell data is passed (impossible to fit kurtosis), 
% feBuildDictionaries will zero out any submitted kurtosis parameters

fe = feBuildDictionaries(fe,Nphi,Ntheta);

%% Encode Connectome in multidimensional tensor
fe = feConnectomeEncoding(fe);
fprintf(['\n[%s] fe-structure Memory Storage:',ByteSize(fe),'\n'],mfilename);

return
