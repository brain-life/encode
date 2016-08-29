function [ Y ] = M_times_w(atoms,voxels,fibers,values,D,w,nTheta,nVoxels)
Y = zeros(nTheta,nVoxels);
for k = 1:length(values)
    %k
    Y(:,voxels(k)) = Y(:,voxels(k)) + D(:,atoms(k))*w(fibers(k))*values(k);
end
Y = Y(:); % vectorization

end