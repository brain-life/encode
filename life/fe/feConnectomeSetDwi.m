function fe = feConnectomeSetDwi(fe,dwiFileName,isrepeat)
% Set all the fields necessary to store the DWI measurements.
%
% This can be used for the dwi measurements use to build the model as well
% as for those used as a second repeat.
%
% Inputs:
%   fe          - an fe structure
%   dwiFileName - the full path to a dwi nifti file
%   isrepeat    - wether the file being set up is the origianl data used to
%                 build LiFE or a repeated measurements used to compute 
%                 statistcs of the quality of fit.        
% Outputs:
%   fe - fe structure with data added in the appropriate fields.
%  
% fe = feConnectomeSetDwi(fe,dwiFileName);   % Set the original data
% fe = feConnectomeSetDwi(fe,dwiFileName,1); % Set the repeated measure
%                                            % data
%
%  Copyright (2020) The Universityof Texas at Austin
%
%  Franco Pestilli frakkopesto@gmail.com and 
%  Cesar F. Caiafa ccaiafa@gmail.com

% Check inputs
if notDefined('isrepeat'), isrepeat=0;end

% Build a tag for the calls to feSet.
if isrepeat, tag = sprintf('repeat');
else,        tag = sprintf('');
end

% Set the ile name in the structure
fe  = feSet(fe, sprintf('dwi%sfile',tag),dwiFileName);

% load the dwi data
dwi = feGet(fe, sprintf('dwi%s',tag));

% Store the bvecs, bvals and the indices for the new file.
fe  = feSet(fe, ...
            sprintf('diffusion bvecs %s',tag), ...
            dwiGet(dwi, 'diffusion bvecs'));
fe  = feSet(fe, ...
            sprintf('diffusion bvals %s',tag), ...
            dwiGet(dwi, 'diffusion bvals'));
fe  = feSet(fe, ...
            sprintf('bvecs indices %s',tag),   ...
            dwiGet(dwi, 'diffusionimagenums'));
            
% Handling multishell data. 
% The original model was only able to handle single-shell. 
% This new version will handle also multishell data.
% 
% To handle multishell data we will store the number of shells encoded
% for the data (1,2,3, etc) and the indices to the bvecs for each shell 
% (this will NOT assume that the number of BVECS is equal for each shell
% (say 30 diffusion directions for each shell).

% We extract the bvalues and find the unique of them and assign them to each shell.
bvals = feGet(fe,'bvals');
[nshells, ~, shellindex] = unique(round(bvals));
fe  = feSet(fe, 'nshells', nshells);
fe  = feSet(fe, 'shellindex', shellindex);

dim = dwi.nifti.dim;
coords = fe.roi.coords;
coords(coords<1)=1;
coords(coords(:,1)>dim(1), 1)=dim(1);
coords(coords(:,2)>dim(2), 2)=dim(2);
coords(coords(:,3)>dim(3), 3)=dim(3);
fe.roi.coords = coords;

% Extract the dwi signal at the coordinates of the connectome
fe  = feSet(fe, ...
            sprintf('diffusion signal image %s',tag), ...
            dwiGet(dwi, 'diffusion signal image',feGet(fe,'roi coords')) );
          
% Extract the non-diffusion direction signal at the coordinates of the
% conenctome
fe  = feSet(fe, sprintf('S0 image %s',tag), ...
            dwiGet(dwi, 'S0 image',feGet(fe,'roi coords')));

if length(isrepeat)>1, keyboard;end
          
% Here I set the dimensions of the dwi file so that I have that
% available everytime when creating a map of parameters.
fe = feSet(fe,sprintf('img size %s',tag),dwiGet(dwi,'size'));

end
