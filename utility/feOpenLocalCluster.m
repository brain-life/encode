function feOpenLocalCluster
% Initializes a local matlab cluster if the parallel matlab toolbox is
% available.
%
%  feOpenLocalCluster
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Initialize a local matlab cluster to speed up some of the processes.
if exist('parpool','file') % check availability of toolbox
   try 
     if isempty(gcp('nocreate')) % if pool not already created
%         if (exist('parcluster','file') == 2)
%            c = parcluster;
%            c.NumWorkers = 12;
%            t = tempname;
%            OK = mkdir(t);
%            if OK
%               c.JobStorageLocation = t;
%            end
%            matlabpool(c);
%         else
%            matlabpool open;     
%         end
        parpool
     else
        % disp('[feOpenLocalCluster] Found Matlab parallel cluster open, not intializing.')
     end
   catch ME
     fprintf('\n[feOpenLocalCluster] Problem intializing the cluster: \n\n %s.\n', ME.message)
     disp('[feOpenLocalCluster] Many computations will be substantially slower, without the parallel toolbox.')
   end
else
      disp('[feOpenLocalCluster] Cannot find the Matlab parallel toolbox, not intializing a cluster.')
      disp('[feOpenLocalCluster] Many computations will be substantially slower, without the parallel toolbox.')
end

end
