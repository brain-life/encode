function [fh, fe] = demo_connectome_encoding()
%% Encode a connectome in multidimensional array (also called tensor).
%
% This demo illustrates how to take as input a tractogrpahy file and a
% diffusion-weighted imaging file and encode them into a multidimensional
% model.
%
%
%  Copyright (2016), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa
%  (CONICET) email: frakkopesto@gmail.com and ccaiafa@gmail.com

% Check if vistasoft is visible on the matlab path.
check = which('vistaRootPath');
if isempty(check)
    error('Vistasoft package not installed. \n Please, download it from https://github.com/vistalab/vistasoft and add it to the Matlab Path');
end

%% Build the file names for the diffusion data, the anatomical MRI.
t1File        = fullfile(feDemoDataPath('HCP3T','sub-105115','anatomy'),  'T1w_acpc_dc_restore_1p25.nii.gz');

%% (1) Load a connectome from disk. 
%
% Here after we refer to connectome to a set of fascicles (streamlines)
% estimated using a tractogrpahy algorithm. This software can read several
% formats for the streamlines:
%  - *.mat and *.pdb from vistasoft; 
%  - *.tck from mrtrix
%
% A connectome generally spans the whole white matter of a brain.
% Streamlines are sets of x,y,z coordinates of 'nodes' whose coordinates
% should be stored in AC-PC (Anterior/Posteriro commissure, RAS,
% Right-Anteriro-Superior format).
%
% We have saved a conenctoem in the demodatsset (URL).
% Here we assume that the dataset has been download and saved in the
% appropriate fodler. See URL for more detaisl on setting up data and file
% paths.

% First, we identify a connectome (a whole brain fascicle group) file on
% disk. In our case this fascicle group was generated using the mrtrix
% toolbox (URL) and the the diffusion-weighted data from the HCP3T data
% set.
fgFileName = fullfile(feDemoDataPath('HCP3T','sub-105115','tractography'), ...
             'dwi_data_b2000_aligned_trilin_csd_lmax10_wm_SD_PROB-NUM01-500000.tck');

%% (2) Identify a DWI file from disk.
%
% Connectomes are the results of a tractography method. For this reason
% they are generated using a diffusion-weighted imaging data set. Diffusion
% weighted images are geenrally saved as either DICOM or NIFTI files. This
% software is compatible with NIFTI files. It can read and write NIFTI
% files using routines from the vistasoft repository.
% 
% Below is the file we used to generate the connectome loaded in (1). We
% will use this DWI data file and encode the DWI data with the connectome. 
dwiFile       = fullfile(feDemoDataPath('HCP3T','sub-105115','dwi'),'dwi_data_b2000_aligned_trilin.nii.gz');

% The final connectome and data astructure will be saved with this name:
feFileName    = 'LiFE_build_model_demo_HCP3T_105115_CSD_PROB';

%% (3) Encode connectome and data in a multidimensional tensor. 
% 
% The multidimensional encoding method that we introduce here has two
% functions: 
% 
% (A) It encodes a fascicles set into a multidimensional array (a
% so-called tensor, to be distinguished from the tensor model geenrally
% used to model the measured diffusion signal).
% 
% (B) It encodes the diffusion-weighted data used to generate the fascicles
% in the connectome into a two-dimensional array, matrix.

% Discretization parameter (this defines the number of orientations encoded
% in Phi, more specifically the size of Phi in mode 1)
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

%% The Phi tensor encodes fascicles' nodes orientation in mode 1 (see Caiafa Pestilli Figure 1)
Number_of_Orientations = feGet(fe,'n atoms');

%% The Phi tensor encodes spatial location of nodes (voxel indices) in mode 2 (see Caiafa Pestilli Figure 1).
Number_of_voxels = feGet(fe,'n voxels');

%% The Phi tensor encodes fascicles identify in mode 3 (see Caiafa Pestilli Figure 1).
Number_of_Fascicles = feGet(fe,'nfibers');

disp(['The size of the sparse tensor Phi is (Na,Nv,Nf) = (',num2str(Number_of_Orientations),',',num2str(Number_of_voxels),',',num2str(Number_of_Fascicles),')'])

%% The precomputed (demeaned) diffusion signals are stored in a Dictionary matrix D. 
% Each column (atom) in the Dictionary corresponds to one spatial orientation of a fascicle's node
%
% To extract the Dictionary matrix you can use feGet.m
D = feGet(fe,'Dictionary');

Number_of_gradient_directions = feGet(fe,'nbvecs');
disp(['The size of the dictionary D is (Ntheta,Na) = (',num2str(Number_of_gradient_directions),',',num2str(Number_of_Orientations),')'])

%% (4) Example of operating on different modes of the tensor:
% In this example we show how to efficiently find fascicles (3rd mode) having a particular
% orientation (1st mode) in the connectome in a neighborhood of a voxel (2nd
% mode).

% For example, in the following we explain how to identify fascicles going
% paralell with axis-z in a particular voxel vecinity.

% Using the function feGetAtoms() we can obtain indices of atoms (columns
% of D), whose orientation is +/-offset degrees appart from the (0,0,1) unit
% vector
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
disp([num2str(length(fascicles_indices)),' fascicles having the orientation ',num2str(main_orient),' in their trajectories, were found'])

% Finally, we generate a visualization of the fascicles and voxels
Visualize_fascicles(fe,fascicles_indices,voxel_indices,'Subset of fascicles meeting orientation criterion')



end

function [] = Visualize_fascicles(fe,fascicles_ind,voxel_ind, fig_name)
colors     = {[.1 .25 .65]};
viewCoords = [0,0];

fg{1}          = feGet(fe,'fibers img');
fg{1}.fibers   = fg{1}.fibers(fascicles_ind);

% plot fascicles
[fh, ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1]);
% Plot voxels
set(gcf,'Color',[1 1 1])
hold on
scatter3(fe.roi.coords(voxel_ind,1)-2,fe.roi.coords(voxel_ind,2)-2,fe.roi.coords(voxel_ind,3)-2,'r')

end


% Local functions to plot the tracts
function [fig_h, light_h] = plotFasciclesNoAnat(fascicles, color, viewCoords, fig_name,tracts_to_clean)
fig_h = figure('name',fig_name,'color','k');
hold on
%set(gca,'visible','off','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78],'Color','w')
set(gca,'visible','off','Color','w')
for iFas  = 1:length(tracts_to_clean)
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
%set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
set(gcf,'Color',[1 1 1])
drawnow


end

