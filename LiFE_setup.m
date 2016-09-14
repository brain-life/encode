root = pwd;
try 
    cd('mexfiles')
    mex -largeArrayDims Mtransp_times_b_mex.c Mtransp_times_b_sub.c -output Mtransp_times_b_mex -DNDEBUG
    
    fprintf('Successfully compiled Mtransp_times_b.\n');
    cd(root)
catch lasterr
    cd(root)
    fprintf('Could not compile Mtransp_times_b.');
    fprintf('WARNING: You can still use the ".m" version but it is extremelly slow!!!!!. You should review how to install an appropriate compiler and use it to generate the mex version on you system. See for example: http://www.mathworks.com/support/compilers/R2015a/index.html ');
    rethrow(lasterr);
end