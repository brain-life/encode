function [FA_profile, SuperFiber, fgResampled] = ...
    dtiFiberGroup_FA_Average(fg, dt, numberOfNodes)
% Average FA across the fibers, along the bundle length
%
% [myFAVals,SuperFiber, fgResampled] = ...
%   dtiFiberGroup_FA_Average(fg, dt, numberOfNodes)
%
% Important assumptions: the fiber bundle is compact. All fibers begin in
% one ROI and end in another. An example of input bundle is smth that
% emerges from clustering or from manualy picking fibers + ends had been
% clipped to ROIs using dtiClipFiberGroupToROIs or smth like that .
%
% OUTPUTS:
%            myFAVals - array with resulting tensor stats. Rows are nodes,
%                     columns correspond to valNames.
%       Superfiber  - a structure describing the mean (core) trajectory and
%                     the spatial dispersion in the coordinates of fibers
%                     within a fiber group.
% fgResampled       - The fiber group that has been resampled to
%                     numberOfNodes and each fiber has been reoriented to
%                     start and end in a consitent location
% Adapted from  dtiFiberGroupPropertyWeightedAverage() in VISTALAB
% Cesar Caiafa (2017)

if notDefined('fg'), error('Fiber group required'); end

% Number of fibers in the whole fiber group
numfibers = size(fg.fibers, 1);

% Should I add checks if fibers need to be reoriented?
% No, dtiComputeSuperFiberRepresentation takes care of that.
% This function will resample the fibers to numberOfNodes and will also
% reorient some fibers, so the notion of "first" and "last" may end up
% converted. [fg] returned in the line below will be resampled and reoriented.  
[SuperFiber, fgResampled] = dtiComputeSuperFiberRepresentation(fg, [], numberOfNodes);

% Compute gradient of SuperFiber
grad = gradient(SuperFiber.fibers{1});
normgrad = sqrt(sum(grad.^2,1));
grad = grad./repmat(normgrad,3,1);

FA_profile = zeros(numberOfNodes,numfibers);

% Estimate node separation
dif = fg.fibers{1}(:,1:end-1)-fg.fibers{1}(:,2:end);
node_separation = mean(sqrt(sum(dif.^2,1)));

parfor node = 1:numberOfNodes
    disp(['node ',num2str(node),'/',num2str(numberOfNodes)])
    xn = SuperFiber.fibers{1}(:,node);
    for fiber = 1:numfibers
        y = fg.fibers{fiber};
        y = y - repmat(xn,1,size(y,2));
        
        %dplane = abs(grad(:,node)'*y); %closest to the perpendicular plane
        d = sqrt(sum(y.^2,1)); % closest to the node
        
        [~, ind] = min(d);
        dplane = abs(grad(:,node)'*y(:,ind)); %closest to the perpendicular plane
        
        if dplane < 10*node_separation
            FA_profile(node,fiber) = GetValFromImage(dt.data, fg.fibers{fiber}(:,ind)', []); %dt.qto_ijk
        else
            FA_profile(node,fiber) = NaN;
        end 
    end
end

FA_profile(FA_profile <= 0) = NaN;

FA_profile = FA_profile';

%FAmean = nanmean(FA_profile,2);
%FAstd = nanstd(FA_profile,0,2);

return
