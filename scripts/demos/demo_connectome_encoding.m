function [fh, fe] = demo_connectome_encoding()
%% Encode a connectome in multidimensional array (also called tensor).
%
% This demo illustrates how to take a tractography file (a full-set of
% streamlines, also called 'fascicles') and associated diffusion-weighted
% imaging data (a NIFTI file plus BVEC/BVAL files used to generate the
% streamlines) and encode them alltogether into a multidimensional tensor
% framework.
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa
%  (CONICET) email: frakkopesto@gmail.com and ccaiafa@gmail.com

%% (0) Check matlab dependencies and path settings.
if ~exist('vistaRootPath.m','file');
    disp('Vistasoft package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/vistalab/vistasoft');
end

if ~exist('mbaComputeFibersOutliers','file')
    disp('ERROR: mba package either not installed or not on matlab path.')
    error('Please, download it from https://github.com/francopestilli/mba')
end

if ~exist('feDemoDataPath.m','file');
    disp('ERROR: demo dataset either not installed or not on matlab path.')
    error('Please, download it from http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz ')
end

%% Build the file names for the diffusion data, the anatomical MRI.
t1File        = fullfile(feDemoDataPath('STN','sub-FP','anatomy'),  ...
                                        't1.nii.gz');

%% (1) Load brain connectome from disk. 
%
% Here after we refer to connectome to a set of fascicles (streamlines)
% estimated using a tractogrpahy algorithm. This software can read several
% formats for the streamlines:
%  - *.mat and *.pdb from vistasoft; 
%  - *.tck from mrtrix
%
% In our framework connectome generally spans the whole white matter of a
% brain and saved. Streamlines are sets of x,y,z coordinates of 'nodes'
% whose coordinates should be stored in AC-PC (Anterior/Posterior
% commissure, RAS, Right-Anteriro-Superior format).
%
% We have saved a connectome in the demodata set (URL:XXXX). Hereafter we
% assume that the dataset has been download and saved in the appropriate
% folder. See URL:XXX for more details on setting up data and file paths.

% First, we identify a connectome (a whole brain streamiles/fascicles
% group) file on disk. In our case this fascicle group was generated using
% the mrtrix toolbox (URL:XXXX) and the the diffusion-weighted data from the
% HCP3T data set (URL:XXXX).
fgFileName = fullfile(feDemoDataPath('STN','sub-FP','tractography'), ...
             'run01_fliprot_aligned_trilin_csd_lmax10_wm_SD_PROB-NUM01-500000.tck');

%% (2) Identify a DWI file from disk.
%
% Sructural brain Connectomes are generated using a tractography method and
% a diffusion-weighted imaging data set. Diffusion weighted images (DWI)
% are generally stored as either DICOM or NIFTI files. This software is
% compatible with NIFTI files. It can read and write NIFTI files using
% routines from the vistasoft repository. The DWI files are expected to be
% preprocessed, this means, at least motion compensated, AC-PC aligned and
% and with artifacts removed.
% 
% Below is the DWI file we used to generate the connectome loaded in (1). We
% will use this DWI data file and encode the DWI data with the connectome. 
dwiFile       = fullfile(feDemoDataPath('STN','sub-FP','dwi'),'run01_fliprot_aligned_trilin.nii.gz');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'LiFE_build_model_demo_STN_FP_CSD_PROB';

%% (3) Encode connectome and data in a multidimensional tensor. 
% 
% The multidimensional encoding method that we introduce here has two
% functions: 
% 
% (A) It encodes a fascicles set (connectome) into a multidimensional array
% (a so-called tensor (to be distinguished from the tensor model geenrally
% used to model the measured diffusion signal).
% 
% (B) It encodes the diffusion-weighted data used to generate the fascicles
% in the connectome into a two-dimensional array, matrix.

% Discretization parameter (this parameter defines the resolution of the
% Phi tensor in describing the orientation of the fascicles in the
% connectome (number of orientations encoded in Phi, more specifically the
% size of Phi in mode 1).
L = 360; 

% The function feConnectomeInit.m collects all the information necessary to
% encode a connectome and build a decomposed LiFE model (see Pestilli et
% al., 2015 and Caiafa & Pestilli Under Review).
%
% Here after we use the function to encode the model first. Below we
% extract the tensor encoding the connectome (Phi).
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFile,t1File,L,[1,0]);

%% The encoding model is comprised by a large, sparse Phi tensor containing the connectome.
%
% To extract the Phi tensor you can use feGet.m
Phi = feGet(fe, 'Phi');

% The Phi tensor encodes fascicles' nodes orientation in mode 1 (see Caiafa
% Pestilli Figure 1)
Number_of_Orientations = feGet(fe,'n atoms');

% The Phi tensor encodes spatial location of nodes (voxel indices) in mode
% 2 (see Caiafa Pestilli Figure 1).
Number_of_voxels = feGet(fe,'n voxels');

% The Phi tensor encodes fascicles identify in mode 3 (see Caiafa Pestilli
% Figure 1).
Number_of_Fascicles = feGet(fe,'nfibers');

disp(['The size of the sparse tensor Phi is (Na,Nv,Nf) = (',num2str(Number_of_Orientations), ...
      ',',num2str(Number_of_voxels),',',num2str(Number_of_Fascicles),')'])

% The precomputed (demeaned) diffusion signals are stored in a Dictionary
% matrix D. Each column (atom) in the Dictionary corresponds to one spatial
% orientation of a fascicle's node
%
% To extract the Dictionary matrix you can use feGet.m
D = feGet(fe,'Dictionary');

Number_of_gradient_directions = feGet(fe,'nbvecs');
disp(['The size of the dictionary D is (Ntheta,Na) = (', ...
          num2str(Number_of_gradient_directions), ...
      ',',num2str(Number_of_Orientations),')'])

%% (4) Example of operating on different modes of the tensor:
% In this example we show how to efficiently find fascicles (3rd mode)
% having a particular orientation (1st mode) in the connectome in a
% neighborhood of a voxel (2nd mode).

% For example, in the following we explain how to identify fascicles going
% paralell with axis-z in a particular voxel vecinity.

% Using the function feGetAtoms() we can obtain indices of atoms (columns
% of D), whose orientation is +/-offset degrees appart from the (0,0,1)
% unit vector
main_orient = [0,0,1]; % Main orientation
offset = 5; % Tolerance in degrees.
atoms_indices = feGetAtoms(fe,main_orient,offset);

% Using the function feGetVoxels() we can obtain indices of voxels in the
% neighborhood of a spatial position ([x,y,z] coordinates)
center_voxel = [76,78,40]; % [x,y,z] coordinates of a center voxel
vicinity_size = 3; 
voxel_indices = feGetVoxels(fe,center_voxel,vicinity_size);

% We restrict our sparse tensor to the orientation (1st mode) meeting the
% criterion (keeping a subset of horizontal slices) for a particular voxel
% vecinity.
Phi_subtensor = Phi(atoms_indices,voxel_indices,:);

% We search for fascicles (3rd mode) having nodes meeting the orientation critierion
% First, we extract the indices of nonzero entries within the subtensor
[inds, ~] = find(Phi_subtensor); % find nonzero entries of subtensor

% Second, we identify fascicle indices for those nonzero entries
fascicles_indices = unique(inds(:,3));
disp([num2str(length(fascicles_indices)),' fascicles having the orientation ', ...
      num2str(main_orient),' in their trajectories, were found'])

% Finally, we generate a visualization of the fascicles and voxels
Visualize_fascicles(fe,fascicles_indices,voxel_indices, ...
                    'Subset of fascicles meeting orientation criterion')

end

% Below are a set of local matlab functions that are sued in this script.
function [] = Visualize_fascicles(fe,fascicles_ind,voxel_ind, fig_name)
% 
% This function is used to visualize the anatomy of part of connectome
% fascicles and a region of interest (ROI). 
% 
% It calls functions from github.com/francopestilli/mba

colors     = {[.1 .25 .65]};
viewCoords = [0,0];

fg{1}          = feGet(fe,'fibers img');
fg{1}.fibers   = fg{1}.fibers(fascicles_ind);

% plot fascicles
[~, ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1]);

% Plot region of interest (ROI), anatomy voxels.
set(gcf,'Color',[1 1 1])
hold on
offset = 1.5;
plot3(fe.roi.coords(voxel_ind,1)-offset, ...
         fe.roi.coords(voxel_ind,2)-offset, ...
         fe.roi.coords(voxel_ind,3)-offset,'ro', ...
         'markerfacecolor','r', 'markersize',15)

end

function [fig_h, light_h] = plotFasciclesNoAnat(fascicles, color, viewCoords, fig_name,tracts_to_clean)
% 
% This function is used to visualize the anatomy of part of connectome
% fascicles. 
%
% It calls functions from github.com/francopestilli/mba

fig_h = figure('name',fig_name,'color','k');
hold on
set(gca,'visible','off','Color','w')
for iFas  = 1:length(tracts_to_clean)
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
set(gcf,'Color',[1 1 1])
drawnow

end