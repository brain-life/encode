function [ fib, vl, vx ] = feConnectionVLSOE( fe, roi1, roi2 )
%% Identify the indices of fibers intersecting ROI1 and ROI1 through the tensor
%
% Inputs:
%     fe   = a fit LiFE fe structure
%     roi1 = a fullpath string to a NIfTI ROI
%     roi2 = a fullpath string to a second NIfTI ROI
%
% Outputs:
%     outFibers = indices of fibers in tensor that intersect both ROIs
%     vx        = internal working data structure, for debugging 
%
% Copyright (2015-2016), Brent McPherson, Franco Pestilli - University of Indiana, Bloomington
%

%% load roi fibers
roi1 = dtiImportRoiFromNifti(roi1);
roi2 = dtiImportRoiFromNifti(roi2);

%% create intermediary data object to manipulate ROI coordinates

% create min / max matrix of coordinates in roi and fe structure
vx.roi1.coords = roi1.coords;
vx.roi2.coords = roi2.coords;

% pull the fe coordinates
vx.fe.coords = fe.roi.coords;
          
% transform roi coordinates into image space for tensor
vx.wrk1.coords = mrAnatXformCoords(fe.life.xform.acpc2img, vx.roi1.coords); 
vx.wrk2.coords = mrAnatXformCoords(fe.life.xform.acpc2img, vx.roi2.coords); 

% round xform coords like dtiExportRoiToNifti
vx.wrk1.coords = ceil(vx.wrk1.coords);
vx.wrk2.coords = ceil(vx.wrk2.coords);

% find the transformed ROI coordinates in the fe object and return the
% indices. These are the fe coords / indices to find fibers in tensor
vx.wrk1.index = ismember(vx.fe.coords, vx.wrk1.coords, 'rows');
vx.wrk2.index = ismember(vx.fe.coords, vx.wrk2.coords, 'rows');

% in order to quickly parse tensor, convert logical to indices
vx.out1.index = find(vx.wrk1.index);
vx.out2.index = find(vx.wrk2.index);

% create sub-index of sptensor
[ inds1, ~ ] = find(fe.life.M.Phi(:, vx.out1.index, :));
[ inds2, ~ ] = find(fe.life.M.Phi(:, vx.out2.index, :));

% pull unique fibers
if size(inds1) > 0
    vx.out1.fibers = unique(inds1(:, 3));
else
    vx.out1.fibers = 0;
end

if size(inds2) > 0
    vx.out2.fibers = unique(inds2(:, 3));
else
    vx.out2.fibers = 0;
end

% find common fibers between ROIs
vx.out.fibers = intersect(vx.out1.fibers, vx.out2.fibers);

%% primary return object

% Fiber Number Weighted Network
fib.fnwn = size(vx.out.fibers, 1);
% the number of fibers between ROIs

% Fiber Density Weighted Network
szROI1 = size(vx.roi1.coords, 1);
szROI2 = size(vx.roi2.coords, 1);
fib.fdwn = (2 * fib.fnwn) / (szROI1 + szROI2);
% weighted fiber counts from: Buchanan et al 2014

% Fiber Length Weighted Network
lengths = cellfun('size', fe.fg.fibers(vx.out.fibers), 2);
fib.flwn = (1 / fib.fnwn) * sum(lengths);
fib.flwn(isnan(fib.flwn)) = 0;
% number of fibers scaled by the summed length

% Fiber Density Corrected by the Fiber Length
fib.fdfl = (2 / (szROI1 + szROI2)) * sum(1 / lengths);
% Hagmann et al, 2008

% Fiber Weight Filtered Count (count of non-zero weighted fibers)
weights = fe.life.fit.weights(vx.out.fibers);
fib.fwfc = sum(weights > 0);
% Pestilli et al, 2014

% Fiber Weight Weighted Network
fib.fwwn = (1 / fib.fnwn) * sum(weights);
fib.fwwn(isnan(fib.fwwn)) = 0;
% Qi et al, 2016

% Fiber Contribution Weighted Network
fib.fcwn = sum(weights);
% Qi et al, 2016

% Fiber Earth Movers Distance / D-Prime - from feVirtual Lesion

% try to run function
try
    [ vl.rmse_wVL, vl.rmse_woVL, vl.nFib_tract, vl.nFib_PN, vl.nVoxels ] = feComputeVirtualLesion(fe, vx.out.fibers);
    fevl = feComputeEvidence(vl.rmse_woVL, vl.rmse_wVL);

    % otherwise, catch output
    fib.femd = fevl.em.mean;
    fib.fsoe = fevl.s.mean;
    
catch
    % if it fails...
    fib.femd = 0;
    fib.fsoe = 0;
    
    % also needs to not be empty
    vl.rmse_wVL = 0;
    vl.rmse_woVL = 0;
    vl.nFib_tract = 0;
    vl.nFib_PN = 0;
    vl.nVoxels = 0;
    
end
% Pestilli et al, 2014

end

