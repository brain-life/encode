function fe = feConnectomeStatistics(fe)
% Compute Curvature and Torsion of fibers
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%

if notDefined('fe'),  error('LiFE (fe = feCreate) struct needed'); end
if ~isfield(fe,'life')
  error('LiFE - the field ''life'' is necessary in the fe structure.')
end

fprintf('\n[%s] Computing Connectome statistics (curvature, torsion, etc) ... ',mfilename); 
tic

nFibers      = feGet(fe,'n fibers');
nTotalNodes = fefgGet(fe.fg,'n total nodes');
fibers = fe.fg.fibers;
imgsize = feGet(fe,'volumesize');
imgsize = imgsize(1:3); % 4th dimension is discarded
nTotalVoxels = prod(imgsize); % including voxels not in the ROI

% Compute fiber statistics
[fiber_id, cur, tor] = fiberStatistics(fibers,nTotalNodes); % this function gives a (1xnTotalNodes) containing the fiber number for each node in tubes

%cur = ones(size(cur));

% Compute voxels
fibers = cell2mat(fibers(:)'); 
voxel_coord = ceil(fibers) + 1;
voxel_id = sub2ind(imgsize, voxel_coord(1,:)', voxel_coord(2,:)', voxel_coord(3,:)');

% Construct Curvature matrix
Curvature = sparse(voxel_id,fiber_id,cur,nTotalVoxels,nFibers);

% Construct Curvature matrix
Torsion = sparse(voxel_id,fiber_id,tor,nTotalVoxels,nFibers);

% The following sparse matrix (Nvoxels x Nfibers) counts the number of
% nodes per fiber per voxel.
A = sparse(voxel_id,fiber_id,ones(nTotalNodes,1),nTotalVoxels,nFibers);

% Construct Indication matrix
unique_ind = unique([voxel_id, fiber_id],'rows');
Indication = sparse(unique_ind(:,1), unique_ind(:,2), ones(size(unique_ind,1),1));

roi_coords = feGet(fe,'roicoords');
roi_ind = sub2ind(imgsize,roi_coords(:,1)',roi_coords(:,2)',roi_coords(:,3)');

% reduce matrices 1st dimension to roi voxels only
Curvature = Curvature(roi_ind,:);
Torsion = Torsion(roi_ind,:);
A = A(roi_ind,:);
Indication = Indication(roi_ind,:);

[voxel_id,fiber_id,vals] = find(Curvature);
nVoxels = length(roi_ind);
a = A(:);
vals = vals./a(sub2ind([nVoxels,nFibers],voxel_id,fiber_id));
Curvature = sparse(voxel_id,fiber_id,vals,nVoxels,nFibers);

[voxel_id,fiber_id,vals] = find(Torsion);
vals = vals./a(sub2ind([nVoxels,nFibers],voxel_id,fiber_id));
Torsion = sparse(voxel_id,fiber_id,vals,nVoxels,nFibers);

fe = feSet(fe,'Curvature',{Curvature,Indication});
fe = feSet(fe,'Torsion',{Torsion,Indication});

clear 'A' 'voxel_id' 'fiber_id' 'vals'

fprintf('took: %2.3fs.\n',toc)

return
end


function [fiber_id, cur, tor] = fiberStatistics(fibers,nTotalNodes)
    fiber_id = zeros(nTotalNodes,1);
    cur = zeros(nTotalNodes,1);
    tor = zeros(nTotalNodes,1);
    node = 1;
    for f = 1:size(fibers,1)
        NumberOfNodes = size(fibers{f},2);
        fiber_id(node:node+NumberOfNodes-1) = f;
        
        [T,N,B,k,t] = mbaFiberProperties(fibers{f}(1,:)', fibers{f}(2,:)', fibers{f}(3,:)');
        
        cur(node:node+NumberOfNodes-1) = k;
        tor(node:node+NumberOfNodes-1) = t;
        
        node = node + NumberOfNodes;
    end
end

function [T,N,B,k,t] = mbaFiberProperties(x,y,z)
% function [T,N,B,k,t] = mbaFiberProperties(x,y,z)
% 
% Computes the Frenet-Serret vectors invariant on the fiber orientation a each node.
%
% This function uses the Frenet-Serret formulas: https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas
%   
% INPUTS: 
%    x, y (and z) are vectors of coordinates (nodes) in the brain. 
%    These vectors represent the x,y and z coordinates of a fibe in #d space. 
%    z can be omitted if the fiber is assuemnd 2D.
%   
%  OUTPUTS:
%  - T: Tangent to the fiber at each node:
%    _
%    T = r' / |r'|
%
%  - N: Normal to the fiber at each node:
%    _ 
%    N =  t' / |T'|
%
%  - B: Binormal to the fiber:
%    _   _   _
%    B = T x N 
%
%  - k: Curvature to the fiber:
%    k = |T'| 
%
%  - t: Torsion of the fiber at each node:
%    t = dot(-B',N)
% 
%    Example:
%    angle = 2*pi*linspace(0,2,200);
%    x = cos(angle);
%    y = sin(a).*sin(a);
%    z = angle/(2*pi);
%    [T,N,B,k,t] = mbaFiberProperties(x,y,z);
%    line(x,y,z), hold on
%    quiver3(x',y',z',T(:,1),T(:,2),T(:,3),'color','r')
%    quiver3(x',y',z',N(:,1),N(:,2),N(:,3),'color','g')
%    quiver3(x',y',z',B(:,1),B(:,2),B(:,3),'color','b')
%    legend('Curve','Tangent','Normal','Binormal')
%    view(0,0)
% 
% FP

if nargin == 2,  z = zeros(size(x)); end

% Make sure input coordinates are oorgnaized as column vectors, we expect columns
x = x(:); y = y(:); z = z(:);

% Compute the gradient of the curve at each node (this can be thoguth fo as the speed of the fiber)
dx = gradient(x);
dy = gradient(y);
dz = gradient(z);
dr = [dx dy dz];

ddx = gradient(dx);
ddy = gradient(dy);
ddz = gradient(dz);
ddr = [ddx ddy ddz];

% Compute the tangent of the fiber at each node
T = dr./mag(dr,3);

% Compute the derivative of the tanget at each node.
dTx =  gradient(T(:,1));
dTy =  gradient(T(:,2));
dTz =  gradient(T(:,3));
dT = [dTx dTy dTz];

% Compute the normal to the fiber trajectory at each node
N = dT./mag(dT,3);

% Compute the binormal to the fiber trajectory at each node
B = cross(T,N);

% Comptue the curvature of the fiber at each node.
% k = mag(dT,1);
k = mag(cross(dr,ddr),1)./((mag(dr,1)).^3);

% Compute the torsion of the fiber at eahc node.
dddx = gradient(ddx); 
dddy = gradient(ddy); 
dddz = gradient(ddz); 
dddr = [dddx dddy dddz];
t = vdot(cross(dr, ddr), dddr) ./ mag(cross(dr, ddr),1).^2;

end

% Helper functions
function N = vdot(A, B) 
% Compute row-wise dot-product of A and B 
 N = zeros(size(A,1),1); 
 for i=1:size(A,1), N(i) = dot(A(i,:), B(i,:)); 
end
end

function N = mag(T,n),
% Compute magnitude of a vector (Nx3)
 N = sum(abs(T).^2,2).^(1/2);
 d = find(N==0); 
 N(d) = eps*ones(size(d));
 N = N(:,ones(n,1));
end

