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
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Check inputs
if notDefined('isrepeat'), isrepeat=0;end

% Build a tag for the calls to feSet.
if isrepeat, tag = sprintf('repeat');
else         tag = sprintf('');
end

% Set the ile name in the structure
fe  = feSet(fe, sprintf('dwi%sfile',tag),dwiFileName);

% load the dwi data
dwi = feGet(fe, sprintf('dwi%s',tag));

% Store the bvecs, bvals and the indices for the new file.
fe  = feSet(fe, sprintf('diffusion bvecs %s',tag), ...
            dwiGet(dwi, 'diffusion bvecs'));
fe  = feSet(fe, sprintf('diffusion bvals %s',tag), ...
            dwiGet(dwi, 'diffusion bvals'));
fe  = feSet(fe, sprintf('bvecs indices %s',tag),   ...
            dwiGet(dwi, 'diffusionimagenums'));

% Extract the dwi signal at the coordinates of the connectome
fe  = feSet(fe, sprintf('diffusion signal image %s',tag), ...
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
