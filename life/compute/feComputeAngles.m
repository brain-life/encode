% This function compute the Bundle-based Minimum Distance (BMD) between two bundles [Garyfallidis et al, 2015]
function [ Angles_matrix] = feComputeAngles(fe, w)
w_ind = find(~isnan(w) & w>0); % indices to fascicles with nnz weight
[Nv] = size(fe.life.M.Phi,2);
[Nf] = size(fe.life.M.Phi,3);

in = [];
jn = [];
ang = [];

parfor v=1:Nv
 %parfor v=1:100000   
     v
    [inds, vals] = find(fe.life.M.Phi(:,v,w_ind));   
    if ~isempty(vals)
        [orientations, unique_fibers] = unify_fibers(inds,fe.life.M.orient);
        Nfv = length(unique_fibers);
        %Nfv = length(vals);
        for n1=1:Nfv
            for n2=n1+1:Nfv
                in = [in, unique_fibers(n1)];
                jn = [jn, unique_fibers(n2)];
                ang = [ang, (180/pi)*acos(dot(orientations(:,n1), orientations(:,n2)))];  
                %in = [in, inds(n1,2)];
                %jn = [jn, inds(n2,2)];
                %ang = [ang, (180/pi)*acos(dot(fe.life.M.orient(:,inds(n1,1)), (fe.life.M.orient(:,inds(n2,1)))))];                 
            end
        end
    end
    
end

Angles_matrix.in=in;
Angles_matrix.jn=jn;
Angles_matrix.ang=ang;
end

% Function that eliminates fibers with zero-weights and repeated fibers assigning a mean
% orientation. It returns a vector with the 3D orientations of fibers in a
% voxel
function [orientations,unique_fibers] = unify_fibers(inds,orient)
unique_fibers = unique(inds(:,2),'rows');
%unique_fibers = intersect(inds(:,2), w_ind);
orientations = zeros(3,size(unique_fibers,1));

for f=1:size(unique_fibers,1)
    [i] = find(inds(:,2)==unique_fibers(f)); % find all the atoms of a same fiber
    orient_sub = orient(:,i);
    orient_sub = mean(orient_sub,2);
    orientations(:,f) = orient_sub./norm(orient_sub);  
end

end
