% This function compute the angles between two tracts. The indices ind1 and
% ind2 indicate the fascicles (3rmode indices) in those tracts,
% respectively
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
function [ Angles] = feComputeAnglesTracts(fe, ind1, ind2)
ind = find(~isnan(fe.life.fit.weights)&fe.life.fit.weights>0); % indices of nnz weights
[Na] = size(fe.life.M.Phi,1); % # of atoms
[Nv] = size(fe.life.M.Phi,2); % # of voxels
[Nf] = size(fe.life.M.Phi,3); % # of fascicles


Phi = fe.life.M.Phi(:,:,ind); % Keeps frontal slices of fascicles with nnz weights only

Phi_tract1 = Phi(:,:,ind1);
[subs, vals] = find(Phi_tract1);
Phi_tract1 = sptensor(subs,ones(size(vals)),size(Phi_tract1)); % put ones in the nnz locations
Phi_tract1 = sparse(subs(:,1),subs(:,2),ones(size(vals)),Na,Nv); % Collapse tensors in 3rd dimension


Phi_tract2 = Phi(:,:,ind2);
[subs, vals] = find(Phi_tract2);
Phi_tract2 = sptensor(subs,ones(size(vals)),size(Phi_tract2)); % put ones in the nnz locations
Phi_tract2 = sparse(subs(:,1),subs(:,2),ones(size(vals)),Na,Nv); % Collapse tensors in 3rd dimension

[o1,v1,vals1] = find(Phi_tract1);
[o2,v2,vals2] = find(Phi_tract2);

Gsubs = [];
Gvals = [];

common_voxels = intersect(v1,v2);

if isempty(common_voxels)
    disp(['Tract do not overlap'])
    Angles=[];
else
    parfor i=1:length(common_voxels)
        disp([num2str(i),'/',num2str(length(common_voxels))]);
        [ind1]=find(v1==common_voxels(i));
        [ind2]=find(v2==common_voxels(i));

        for j1 = 1:length(o1(ind1))
            for j2 = 1:length(o2(ind2))
                 Gsubs = [Gsubs ; [o1(ind1(j1)), o2(ind2(j2)), common_voxels(i)]];
                 Gvals = [Gvals ; [vals1(ind1(j1))*vals2(ind2(j2))]];
            end
        end 
    end

    Gangles = (180/pi)*acos(dot(fe.life.M.orient(:,Gsubs(:,1)),fe.life.M.orient(:,Gsubs(:,2)),1));
    Gangles = real(Gangles');
    
    ind_higher90 = find(Gangles>90);
    Gangles(ind_higher90) = 180*ones(size(ind_higher90)) - Gangles(ind_higher90);
    
    Npoints = sum(Gvals);
    Angles = zeros(Npoints,1);

    p=1;
    for j=1:size(Gvals)
        Rep = Gvals;
        while Rep
            Angles(p) = Gangles(j);
            p = p + 1;
            Rep = Rep - 1;
        end
    end
    
    Angles_sym = Angles;
    Angles_sym = 180*ones(size(Angles)) - Angles_sym;
    
    Angles = [Angles; Angles_sym];
    
end


end
