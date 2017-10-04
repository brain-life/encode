function [val1] = GetValFromImage(Im, coords, xform)
%
% val = GetValFromImage(Im, coords, [xform])
%
% Adapted from dtiGetValFromImage in VISTASOFT (interpolation was eliminated)
% Cesar Caiafa (2017)

if(~exist('xform','var') || isempty(xform))
    xform = eye(4);
end

if(size(coords,2)~=3) coords = coords'; end
if(size(coords,2)~=3) error('coords must be an Nx3 array!'); end

if(~all(all(xform==eye(4))))
    coords = mrAnatXformCoords(xform, coords);
end

[xyz] = round(coords) + ones(size(coords)); % 

val1 = Im(xyz(:,1),xyz(:,2),xyz(:,3));

return
