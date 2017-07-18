function checkMexCompiled(varargin)
%CHECKMEXCOMPILED Check if mex file is compiled for system
% 
%   IOSR.GENERAL.CHECKMEXCOMPILED(SOURCE_FILE) checks whether a mex source
%   file SOURCE_FILE is compiled for the current operating system OR
%   whether the source file has been modified since it was compiled. It is
%   compiled if it does not pass these tests (to the same directory as the
%   source file). SOURCE_FILE must be a string that is the name of a source
%   file on the MATLAB search path.
% 
%   IOSR.GENERAL.CHECKMEXCOMPILED(OPTIONS,...,SOURCE_FILE) passes the
%   script switches in OPTIONS to the mex compiler, one argument per
%   switch.
% 
%   Example
% 
%       % check function compiled, with debugging info, and
%       % with large-array-handling API
%       iosr.general.checkMexCompiled('-g','-largeArrayDims','myfun.c')
% 
%   See also MEX.

%   Copyright 2016 University of Surrey.

    source_file{1} = varargin{end-1};
    source_file{2} = varargin{end};

    % Check input filename
    assert(ischar(source_file{1}),'source_file must be a string')
    assert(ischar(source_file{2}),'source_file must be a string')

    % Check extension is specified
    assert(~isempty(strfind(source_file{1},'.')),'source_file: no file extension specified')
    assert(~isempty(strfind(source_file{2},'.')),'source_file: no file extension specified')

    % Locate source file
    [pathstr{1},name{1},ext{1}] = fileparts(which(source_file{1}));
    [pathstr{2},name{2},ext{2}] = fileparts(which(source_file{2}));

    filename{1} = [pathstr{1} filesep name{1} ext{1}]; % Create filename
    mexfilename = [pathstr{1} filesep name{1} '.' mexext]; % Deduce mex file name based on current platform

    if strcmp(pathstr{1},'') || strcmp(pathstr{2},'')% source file not found
        error([source_file ': not found'])
    elseif exist(mexfilename,'file')~=3 || get_mod_date(mexfilename)<get_mod_date(filename{1})
         % if source file does not exist or it was modified after the mex file
        disp(['Compiling "' name ext '".'])
        d = cd;
        cd(pathstr{1})
        % compile, with options if appropriate
        try
            if length(varargin)>1
                options = varargin{1:end-2};
                mex(options,source_file{1},source_file{2})
            else
                mex(source_file,source_file{1},source_file{2})
            end
             fprintf('Function %s successfuly compiled\n',source_file{1});
        catch lasterr
            fprintf('ERROR: Could not compile function %s. \n',source_file{1});
            error('Please, install an appropriate compiler and run the scipts again (see http://www.mathworks.com/support/compilers). Mex files will be generated automatically. \n');
        end
        
        cd(d)
    end

end

function datenum = get_mod_date(file)
%GET_MOD_DATE get file modified date

    d = dir(file);
    datenum = d.datenum;

end
