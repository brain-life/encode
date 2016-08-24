function [ w ] = Mtransp_times_b_NOloop(atoms,voxels,fibers,values,D,Y)
%  Matlab only version of matrix multiply between factorized LiFE model and diffusion to
%  compute the weights.
%  
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
w  = dot(D(:,atoms),Y(:,voxels)).*values';
w = accumarray(fibers,w');

end



