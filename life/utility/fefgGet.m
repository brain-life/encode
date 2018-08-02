function val = fefgGet(fg,param,varargin)
% Get values from a fiber group structure
%
%  val = fefgGet(fg,param,varargin)
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
%
%------------
% Returns the length of each fiber in the fiber group
% flen = fefgGet(fg,'length')
%------------
% Compute the Westin shape parameters (linearity and planarity)
% for all the fibers in the fiber group.
% Requires a dt structure.
% westin = fefgGet(fg,'westinshape',dt)
% westin = fefgGet(fg,'westinshape',eigenvals)
% westin = fefgGet(fg,'westinshape',dtFileName)
%-------------
% Compute eigenvalues for the fibers in a fiber group.
% It requires a dt structure.
% eigenvals = fefgGet(fg,'eigenvals',dt)
% eigenvals = fefgGet(fg,'eigenvals',dtFileName)
%-------------
% Compute the 6 parameters for the tensors at each node in each fiber.
% It requires a dt structure.
% dt6 = fefgGet(fg,'eigenvals',dt)
% dt6 = fefgGet(fg,'eigenvals',dtFileName)
%-------------
% Compute the radial and axial ADC for all the fibers in the fiber group.
% Requires a dt structure.
% adc = fefgGet(fg,'adc',dt)
% adc = fefgGet(fg,'adc',eigenvals)
% adc = fefgGet(fg,'adc',dtFileName)
%-------------
% Compute the FA for all the fibers in the fiber group.
% Requires a dt structure.
% fa = fefgGet(fg,'fa',dt)
% fa = fefgGet(fg,'fa',eigenvals)
% fa = fefgGet(fg,'fa',dtFileName)
%-------------
% Compute the Mean Diffusivity for all the fibers in the fiber group.
% Requires a dt structure.
% md = fefgGet(fg,'md',dt)
% md = fefgGet(fg,'md',eigenvals)
% md = fefgGet(fg,'md',dtFileName)
%-------------
% Compute the Radial Diffusivity for all the fibers in the fiber group.
% Requires a dt structure.
% rd = fefgGet(fg,'rd',dt)
% rd = fefgGet(fg,'rd',eigenvals)
% rd = fefgGet(fg,'rd',dtFileName)
%-------------
% Compute the Axial Diffusivity for all the fibers in the fiber group.
% Requires a dt structure.
% ad = fefgGet(fg,'ad',dt)
% ad = fefgGet(fg,'ad',eigenvals)
% ad = fefgGet(fg,'ad',dtFileName)
%-------------
% Compute the radial and axial ADC for all the fibers in the fiber group.
% Requires a dt structure.
% adc = fefgGet(fg,'adc',dt)
% adc = fefgGet(fg,'adc',eigenvals)
% adc = fefgGet(fg,'adc',dtFileName)
%-------------
%
% Parameters
% General
%    'name'
%    'type'
%    'colorrgb'
%    'thickness'
%    'visible'
%
% Fiber related
%    'nfibers'- Number of fibers in this group
%    'nodes per fiber'  - Number of nodes per fiber.
%    'fibers' - Fiber coordinates
%    'fibernames'
%    'fiberindex'
% Compute the Mean Diffusivity for all the fibers in the fiber group.
% Requires a dt structure.
% fa = fefgGet(fg,'fa',dt)
% fa = fefgGet(fg,'fa',eigenvals)
%
% ROI and image coord related
%    'unique image coords'
%    'nodes to imagecoords' -
%    'voxel2fiber node pairs' - For each roi coord, an Nx2 matrix of
%        (fiber number,node number)
%    'nodes in voxels' - Nodes inside the voxels of roi coords
%    'voxels in fg'  - Cell array of the roiCoords touched by each fiber
%    'voxels2fibermatrix' - Binary matrix (voxels by fibers).  1s when a
%        fiber is in a voxel of the roiCoords (which are, sadly, implicit).
%
% Tensor and tractography related
%      'tensors'     - Tensors for each node
%
%
% See also: fgCreate; fgSet

val = [];

switch strrep(lower(param),' ','')
  
  % Basic fiber parameters
  case 'name'
    val = fg.name;
  case 'type'  % Should always be fibergroup
    val = fg.type;
    
    % Fiber visualization settings.
  case 'colorrgb'
    val = fg.colorRgb;
  case 'thickness'
    val = fg.thickness;
  case 'visible'
    val = fg.visible;
    
    % Simple fiber properties --
  case {'fibers'}
    % val = fefgGet(fg,'fibers',fList);
    %
    % Returns a 3xN matrix of fiber coordinates corresponding to the
    % fibers specified in the integer vector, fList.  This differs from
    % the dtiH (mrDiffusion) representation, where fiber coordinates
    % are stored as a set of cell arrays for each fiber.
    if ~isempty(varargin)
      list = varargin{1};
      val = cell(length(list),1);
      for ii=1:length(list)
        val{ii} = fg.fibers{ii};
      end
    else
      val = fg.fibers;
    end
  case 'fibernames'
    val = fg.fiberNames;
  case 'fiberindex'
    val = fg.fiberIndex;
  case 'nfibers'
    val = length(fg.fibers);
  case {'nodesperfiber','nsamplesperfiber','nfibersamples'}
    % fefgGet(fg,'n samples per fiber ')
    % How many samples per fiber.  This is about equal to
    % their length in mm, though we need to write the fiber lengths
    % routine to actually calculate this.
    nFibers = fefgGet(fg,'n fibers');
    val = zeros(1,nFibers);
    for ii=1:nFibers
      val(ii) = length(fg.fibers{ii});
    end
    
    % Fiber group (subgroup) properties.
    % These are used when we classify fibers into subgroups.  We should
    % probably clean up this organization which is currently
    %
    %   subgroup - length of fibers, an index of group identity
    %   subgroupNames()
    %    .subgroupIndex - Probably should go away and the index should
    %          just be
    %    .subgroupName  - Probably should be moved up.
    %
  case {'ngroups','nsubgroups'}
    val = length(fg.subgroupNames);
  case {'groupnames'}
    val = cell(1,fefgGet(fg,'n groups'));
    for ii=1:nGroups
      val{ii} = fg.subgroupNames(ii).subgroupName;
    end
    
    % DTI properties
  case 'tensors'
    val = fg.tensors;
    
    % Fiber to coord calculations
  case {'imagecoords'}
    % c = fefgGet(fgAcpc,'image coords',fgList,xForm);
    % c = fefgGet(fgAcpc,'image coords',fgList,xForm);
    %
    % Return the image coordinates of a specified list of fibers
    % Returns a matrix that is fgList by 3 of the image coordinates for
    % each node of each fiber.
    %
    % Fiber coords are represented at fine resolution in ACPC space.
    % These coordinates are rounded and in image space
    if ~isempty(varargin)
      fList = varargin{1};
      if length(varargin) > 1
        xForm = varargin{2};
        % Put the fiber coordinates into image space
        fg = dtiXformFiberCoords(fg,xForm);
      end
    else
      % In this case, the fiber coords should already be in image
      % space.
      nFibers = fefgGet(fg,'n fibers');
      fList = 1:nFibers;
    end
    
    % Pull out the coordinates and round them.  These are in image
    % space.
    nFibers = length(fList);
    val = cell(1,nFibers);
    if nFibers == 1
      val = round(fg.fibers{fList(1)}')+1;

    else
     feOpenLocalCluster
     parfor ii=1:nFibers
         val{ii} = round(fg.fibers{fList(ii)}')+1;
     end
    end

  case {'uniqueimagecoords'}
    %   coords = fefgGet(fgIMG,'unique image coords');
    %
    % The fg input must be in IMG space.
    %
    % Returns the unique image coordinates of all the fibers as an Nx3
    % matrix of integers.
    % val = round(horzcat(fg.fibers{:})'); 
    val = round(horzcat(fg.fibers{:})')+1;
    val = unique(val,'rows');
  
    case {'allimagecoords'}
    %   coords = fefgGet(fgIMG,'all image coords');
    %
    % The fg input must be in IMG space.
    %
    % Returns all image coordinates of all the fibers as an Nx3
    % matrix of integers.
    % val = round(horzcat(fg.fibers{:})'); 
    val = round(horzcat(fg.fibers{:})')+1;
    
  case {'uniqueacpccoords'}
    %   coords = fefgGet(fg,'unique acpc coords');
    %
    % The fg input must be in ACPC space.
    %
    % Returns the unique ACPC coordinates of all the fibers as an Nx3
    % matrix of integers.
    val = fefgGet(fg, 'uniqueimagecoords') ;

  case {'nodes2voxels'}
    %   nodes2voxels = fefgGet(fgImg,'nodes2voxels',roiCoords)
    %
    % The roiCoords are a matrix of Nx3 coordinates.  They describe a
    % region of interest, typically in image space or possibly in acpc
    % space.
    %
    % We return a cell array that is a mapping of fiber nodes to voxels in
    % the roi.  The roi is specified as an Nx3 matrix of coordinates.
    % The returned cell array, nodes2voxels, has the same number of
    % cells as there are fibers.
    %
    % Unlike the fiber group cells, which have a 3D coordinate of each
    % node, this cell array has an integer that indexes the row of
    % roiCoords that contains the node. If a node is not in any of the
    % roiCoords, the entry in node2voxels{ii} for that node is zero.
    % This means that node is outside the 'roiCoords'.
    %
    % Once again: The cell nodes2voxels{ii} specifies whether each
    % node in the iith fiber is inside a voxel in the roiCoords.  The
    % value specifies the row in roiCoords that contains the node.
    %
    if isempty(varargin), error('roiCoords required');
    else
      roiCoords = varargin{1};
    end
    fprintf('[%s] Computing nodes-to-voxels..',mfilename)

    % Find the roiCoord for each node in each fiber.
    nFiber = fefgGet(fg,'n fibers');
    val    = cell(nFiber,1);
    
   parfor ii=1:nFiber
     % if ~mod(ii,200), fprintf('%d ',ii); end
     % Node coordinates in image space
     nodeCoords = fefgGet(fg,'image coords',ii);
     
     % The values in loc are the row of the coords matrix that contains
     % that sample point in a fiber.  For example, if the number 100 is
     % in the 10th position of loc, then the 10th sample point in the
     % fiber passes through the voxel in row 100 of coords.
     [~, val{ii}] = ismember(nodeCoords, roiCoords, 'rows');  % This operation is slow.
    end    
    fprintf(' took: %2.3fminutes.\n',toc/60)

  case {'voxel2fibernodepairs','v2fn'}
    % voxel2FNpairs = fefgGet(fgImg,'voxel 2 fibernode pairs',roiCoords);
    % voxel2FNpairs = fefgGet(fgImg,'voxel 2 fibernode pairs',roiCoords,nodes2voxels);
    %
    % The return is a cell array whose size is the number of voxels.
    % The cell is a Nx2 matrix of the (fiber, node) pairs that pass
    % through it.
    %
    % The value N is the number of nodes in the voxel.  The first
    % column is the fiber number.  The second column reports the indexes
    % of the nodes for each fiber in each voxel.
    fprintf('[%s] Computing fibers/nodes pairing in each voxel..\n',mfilename)
   
    if (length(varargin) < 1), error('Requires the roiCoords.');
    else
      roiCoords = varargin{1};
      nCoords   = size(roiCoords,1);
    end
    if length(varargin) < 2
      % We assume the fg and the ROI coordinates are in the same
      % coordinate frame. The following line gives the list of fibers
      % containing only the nodes that are included in the ROI
      tic,nodes2voxels    = fefgGet(fg,'nodes 2 voxels',roiCoords);
    else nodes2voxels = varargin{2};
    end
    
    nFibers      = fefgGet(fg,'nFibers');
    voxelsInFG   = fefgGet(fg,'voxels in fg',nodes2voxels);      

    tic,roiNodesInFG = fefgGet(fg,'nodes in voxels',nodes2voxels);
    tic, val = cell(1,nCoords);
    for thisFiber=1:nFibers
      voxelsInFiber = voxelsInFG{thisFiber};   % A few voxels, in a list
      nodesInFiber  = roiNodesInFG{thisFiber}; % The corresponding nodes
      
      % Then add a row for each (fiber,node) pairs that pass through
      % the voxels for this fiber.
      for jj=1:length(voxelsInFiber)
        thisVoxel = voxelsInFiber(jj);
        % Print out roi coord and fiber coord to verify match
        % roiCoords(thisVoxel,:)
        % fg.fibers{thisFiber}(:,nodesInFiber(jj))
        % Would horzcat be faster?
        val{thisVoxel} = cat(1,val{thisVoxel},[thisFiber,nodesInFiber(jj)]);
      end
    end
    fprintf('[%s] fiber/node pairing completed in: %2.3fs.\n',mfilename, toc)
  
  case {'voxel2fibers','fiberdensity'}
    % voxel2FNpairs = fefgGet(fgImg,'fiber density',roiCoords);
    %
    % The return is a cell array whose size is the number of voxels.
    % The cell is a Nx1 matrix of fiber counts. How many fibers in each
    % voxel.
    %
    fprintf('[%s] Computing fiber density in each voxel...\n',mfilename)
   
    if (length(varargin) < 1), 
      roiCoords = fefgGet(fg,'uniqueimagecoords');
      fprintf('[%s] Computing white-matter volume ROI from fibers coordinates.\n',mfilename) 
      fprintf('          Assuming fiber coordinates in IMG space.\n')
    else
      roiCoords = varargin{1};
    end
    
    if length(varargin) < 2
      % We assume the fg and the ROI coordinates are in the same
      % coordinate frame.
      tic,nodes2voxels= fefgGet(fg,'nodes 2 voxels',roiCoords);
    else nodes2voxels = varargin{2};
    end
    
    nCoords      = size(roiCoords,1);
    nFibers      = fefgGet(fg,'nFibers');
    voxelsInFG   = fefgGet(fg,'voxels in fg',nodes2voxels);      

    tic, fibersInVox = cell(1,nCoords);
    for thisFiber=1:nFibers
      voxelsInFiber = voxelsInFG{thisFiber};   % A few voxels, in a list
      
      % Then add a row for each (fiber,node) pairs that pass through
      % the voxels for this fiber.
      for jj=1:length(voxelsInFiber)
        thisVoxel = voxelsInFiber(jj);
        % Print out roi coord and fiber coord to verify match
        % roiCoords(thisVoxel,:)
        % fg.fibers{thisFiber}(:,nodesInFiber(jj))
        % Would horzcat be faster?
        fibersInVox{thisVoxel} = cat(1,fibersInVox{thisVoxel},thisFiber);
      end
    end
    
    % Now create the fiber density, the unique fibers in each voxel
    val = zeros(length(fibersInVox),1);
    for ivx = 1:length(fibersInVox)
      val(ivx) = length(unique(fibersInVox{ivx}));
    end
    
    fprintf('[%s] fiber density completed in: %2.3fs.\n',mfilename, toc)
  
  case {'uniquefibersinvox'}
    % uniquefvx = fefgGet(fgImg,'uniquefibrsinvox',roiCoords);
    %
    % The return is a vector whose size is the number of voxels containing
    % the unique fibers in each voxel
    %
    fprintf('[%s] Computing the unique fibers in each voxel...\n',mfilename)
    tic
    if (length(varargin) < 1), error('Requires the roiCoords.');
    else
      roiCoords = varargin{1};
      nCoords   = size(roiCoords,1);
    end
    
    % We assume the fg and the ROI coordinates are in the same
    % coordinate frame.
    nodes2voxels    = fefgGet(fg,'nodes 2 voxels',roiCoords);
    
    nFibers      = fefgGet(fg,'nFibers');
    voxelsInFG   = fefgGet(fg,'voxels in fg',nodes2voxels);      

    fibersInVox = cell(1,nCoords);
    for thisFiber=1:nFibers
      voxelsInFiber = voxelsInFG{thisFiber};   % A few voxels, in a list
      
      % Then add a row for each (fiber,node) pairs that pass through
      % the voxels for this fiber.
      if ~isempty(voxelsInFiber)
          for jj=1:length(voxelsInFiber)
              thisVoxel = voxelsInFiber(jj);
              fibersInVox{thisVoxel} = cat(1,fibersInVox{thisVoxel},thisFiber);
          end
      end
    end
    
    % Now create find the unique fibers in each voxel
    val = cell(length(fibersInVox),1);
    for ivx = 1:length(fibersInVox)
      val{ivx} = (unique(fibersInVox{ivx}));
    end
    fprintf('[%s] done computing unique fibers in each voxel: %2.3fs.\n',mfilename, toc)

  case {'nodesinvoxels'}
    % nodesInVoxels = fefgGet(fg,'nodes in voxels',nodes2voxels);
    %
    % This cell array is a modified form of nodes2voxels (see above).
    % In that cell array every node in every fiber has a number
    % referring to its row in roiCoords, or a 0 when the node is not in
    % any roiCoord voxel.
    %
    % This cell array differs only in that the 0s removed.  This
    % is used to simplify certain calculations.
    %
    if (length(varargin) <1), error('Requires nodes2voxels cell array.'); end
    tic,fprintf('[%s] Computing nodes-in-voxels..',mfilename)
    nodes2voxels = varargin{1};
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
 
    % For each fiber, this is a list of the nodes that pass through
    % a voxel in the roiCoords
    for ii = 1:nFibers
      % For each fiber, this is a list of the nodes that pass through
      % a voxel in the roiCoords
      lst = (nodes2voxels{ii} ~= 0);
      val{ii} = find(lst);
    end
    fprintf(' took: %2.3fs.\n',toc)

  case 'voxelsinfg'
    % voxelsInFG = fefgGet(fgImg,'voxels in fg',nodes2voxels);
    %
    % A cell array length n-fibers. Each cell has a list of the voxels
    % (rows of roiCoords) for a fiber.
    %
    % This routine eliminates the 0's in the nodes2voxels lists.
    %
    if length(varargin) < 1, error('Requires nodes2voxels cell array.'); end
    tic,fprintf('[%s] Computing voxels-in-fg..',mfilename)

    nodes2voxels = varargin{1};
    nFibers = fefgGet(fg,'nFibers');
    val = cell(1,nFibers);
 
    for ii = 1:nFibers
      % These are the nodes that pass through a voxel in the
      % roiCoords
      lst = (nodes2voxels{ii} ~= 0);
      val{ii} = nodes2voxels{ii}(lst);
    end
    fprintf(' took: %2.3fs.\n',toc)

  case {'voxels2fibermatrix','v2fm'}
    %   v2fm = fefgGet(fgImg,'voxels 2 fiber matrix',roiCoords);
    % Or,
    %   v2fnPairs = fefgGet(fgImg,'v2fn',roiCoords);
    %   v2fm = fefgGet(fgImg,'voxels 2 fiber matrix',roiCoords, v2fnPairs);
    %
    % mrvNewGraphWin; imagesc(v2fm)
    %
    % Returns a binary matrix of size Voxels by Fibers.
    % When voxel ii has at least one node from fiber jj, there is a one
    % in v2fm(ii,jj).  Otherwise, the entry is zero.
    %
    
    % Check that the fg is in the image coordspace:
    if isfield(fg, 'coordspace') && ~strcmp(fg.coordspace, 'img')
      error('Fiber group is not in the image coordspace, please xform');
    end
    
    if isempty(varargin), error('roiCoords required');
    else
      roiCoords = varargin{1};
      nCoords   = size(roiCoords,1);
      if length(varargin) < 2
        v2fnPairs = fefgGet(fg,'v2fn',roiCoords);
      else
        v2fnPairs = varargin{2};
      end
    end
    
    % Allocate matrix of voxels by fibers
    val = zeros(nCoords,fefgGet(fg,'n fibers'));
    
    % For each coordinate, find the fibers.  Set those entries to 1.
    for ii=1:nCoords
      if ~isempty(v2fnPairs{ii})
        f = unique(v2fnPairs{ii}(:,1));
      end
      val(ii,f) = 1;
    end
    
  case {'fibersinroi','fginvoxels','fibersinvoxels'}
    % fList = fefgGet(fgImg,'fibersinroi',roiCoords);
    %
    % v2fn = fefgGet(fgImg,'v2fn',roiCoords);
    % fList = fefgGet(fgImg,'fibersinroi',roiCoords,v2fn);
    %
    % Returns an integer vector of the fibers with at least
    % one node in a region of interest.
    %
    % The fg and roiCoords should be in the same coordinate frame.
    %
    if isempty(varargin), error('roiCoords required');
    elseif length(varargin) == 1
      roiCoords = varargin{1};
      v2fnPairs = fefgGet(fg,'v2fn',roiCoords);
    elseif length(varargin) > 1
      roiCoords = varargin{1};
      v2fnPairs = varargin{2};
    end
    
    val = []; nCoords = size(roiCoords,1);
    for ii=1:nCoords
      if ~isempty(v2fnPairs{ii})
        val = cat(1,val,v2fnPairs{ii}(:,1));
      end
    end
    val = sort(unique(val),'ascend');
    
  case {'coordspace','fibercoordinatespace','fcspace'}
    % In some cases, the fg might contain information telling us in which
    % coordinate space its coordinates are set. This information is set
    % as a struct. Each entry in the struct can be either a 4x4 xform
    % matrix from the fiber coordinates to that space (with eye(4) for
    % the space in which the coordinates are defined), or (if the xform
    % is not know) an empty matrix.
    
    cspace_fields = fields(fg.coordspace);
    val = [];
    for f=1:length(cspace_fields)
      this_field = cspace_fields{f};
      if isequal(getfield(fg.coordspace, this_field), eye(4))
        val = this_field;
      end
    end
  
  case {'length'}
    % Returns the length of each fiber in the fiber group
    % flen = fefgGet(fg,'length')
    
    % Measure the step size of the first fiber. They *should* all be the same!
    stepSize = mean(sqrt(sum(diff(fg.fibers{1},1,2).^2)));
    
    % Check that they are
    %stepSize2 = mean(sqrt(sum(diff(fg.fibers{2},1,2).^2)));
    %assertElementsAlmostEqual(stepSize,stepSize2,'relative',.001,.001);
    
    % Estimate the length for the fibers, as well mean, min and max
    val  = cellfun('length',fg.fibers)*stepSize;
  
  case {'dt6'}
    % Compute the 6 parameters for the tensors at each node in each fiber.
    % It requires a dt structure.
    % dt6 = fefgGet(fg,'dt6',dt)  
    % dt6 = fefgGet(fg,'dt6',dtFileName)

    if ~ischar(varargin{1})    
      if isstruct(varargin{1})
        % A dti structure was passed.
        dt = varargin{1};
      else
        error('[%s] A ''dt'' structure is necessary.', mfilename)
      end
    else % The string should be a path to a dt.mat file.
      dt = dtiLoadDt6(varargin{1});
    end
    
    % Get the dt6 tensors for each node in each fiber.
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    xform   = inv(dt.xformToAcpc);
    dt6     = dt.dt6; clear dt
    for ii = 1:nFibers
      % Get the fiber coordinates.
      coords = fg.fibers{ii}'; % Assumed in ACPC
      
      [val1,val2,val3,val4,val5,val6] = ...
        dtiGetValFromTensors(dt6, coords, xform,'dt6','nearest');
      
      % Build a vector of tesnsors' eigenvalues
      val{ii} = [val1,val2,val3,val4,val5,val6];
    end
    
    
  case {'eigenvals'}
    % Compute eigenvalues for the fibers in a fiber group.
    % It requires a dt structure.
    % eigenvals = fefgGet(fg,'eigenvals',dt)
    % eigenvals = fefgGet(fg,'eigenvals',dt6FileName)
    if ~ischar(varargin{1})
      if ( ~isstruct(varargin{1}) )
        dt6 = fefgGet(fg,'dt6',varargin{1});
      elseif ( size(varargin{1},1)==6 )
        dt6 = varargin{1};
      else
        error('[%s] Either a ''dt'' structure or a ''dt6'' vector of tensor coefficients is necessary.\n',mfilename)
      end
    else % The string should be a path to a dt6.mat file.
      dt6 = fefgGet(fg,'dt6',dtiLoadDt6(varargin{1}));
    end
    
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    for ii = 1:nFibers
      
      % We now have the dt6 data from all of the fibers.  We extract the
      % directions into vec and the eigenvalues into val.  The units of val are
      % um^2/sec or um^2/msec .. somebody answer this here, please.
      [~,val{ii}] = dtiEig(dt6{ii});
    end
    
    for  ii = 1:nFibers
      % Some of the ellipsoid fits are wrong and we get negative eigenvalues.
      % These are annoying. If they are just a little less than 0, then clipping
      % to 0 is not an entirely unreasonable thing. Maybe we should check for the
      % magnitude of the error?
      nonPD = find(any(val{ii}<0,2), 1);
      if(~isempty(nonPD))
        %fprintf('\n NOTE: %d fiber points had negative eigenvalues. These will be clipped to 0..\n',length(nonPD));
        val{ii}(val{ii}<0) = 0;
      end
      
      threeZeroVals = (sum(val{ii}, 2) == 0);
      %if ~isempty (threeZeroVals)
      %  fprintf('\n NOTE: %d of these fiber points had all three negative eigenvalues. These will be excluded from analyses\n', length(threeZeroVals));
      %end
      val{ii}(threeZeroVals, :)=[];
    end
    
  case {'fa'}
    % Compute the FA for all the fibers in the fiber group.
    % Requires a dt structure.
    % fa = fefgGet(fg,'fa',dt)
    % fa = fefgGet(fg,'fa',eigenvals)
    % fa = fefgGet(fg,'fa',dtFileName)
    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1})
      % A dti structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    
    for ii = 1:nFibers    
      [val{ii},~,~,~] = dtiComputeFA(eigenvals{ii});
    end
    
  case {'md'}
    % Compute the Mean Diffusivity for all the fibers in the fiber group.
    % Requires a dt structure.
    % md = fefgGet(fg,'md',dt)
    % md = fefgGet(fg,'md',eigenvals)
    % md = fefgGet(fg,'md',dtFileName)

    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1})
      % A dti structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    
    for ii = 1:nFibers    
      [~,val{ii},~,~] = dtiComputeFA(eigenvals{ii});
    end
       
  case {'ad'}
    % Compute the Axial Diffusivity for all the fibers in the fiber group.
    % Requires a dt structure.
    % ad = fefgGet(fg,'ad',dt)
    % ad = fefgGet(fg,'ad',eigenvals)
    % ad = fefgGet(fg,'ad',dtFileName)

    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1}) 
      % A dti structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    
    for ii = 1:nFibers    
      [~,~,~,val{ii}] = dtiComputeFA(eigenvals{ii});
    end
        
  case {'rd'}
    % Compute the Radial Diffusivity for all the fibers in the fiber group.
    % Requires a dt structure.
    % rd = fefgGet(fg,'rd',dt)
    % rd = fefgGet(fg,'rd',eigenvals)
    % rd = fefgGet(fg,'rd',dtFileName)

    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1}) 
      % A dti structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    
    for ii = 1:nFibers    
      [~,~,val{ii},~] = dtiComputeFA(eigenvals{ii});
    end
    
  case {'adc'}
    % Compute the radial and axial ADC for all the fibers in the fiber group.
    % Requires a dt structure.
    % adc = fefgGet(fg,'adc',dt)
    % adc = fefgGet(fg,'adc',eigenvals)
    % adc = fefgGet(fg,'adc',dtFileName)
    
    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1}) 
      % A dti structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);
    for ii = 1:nFibers
      [~,~,val{ii}.radial, val{ii}.axial] = dtiComputeFA(eigenvals{ii});
    end
    
  case {'westinshape'}
    % Compute the Westin shape parameters (linearity and planarity)
    % for all the fibers in the fiber group.
    % Requires a dt structure.
    % westin = fefgGet(fg,'westinshape',dt)
    % westin = fefgGet(fg,'westinshape',eigenvals)
    % westin = fefgGet(fg,'westinshape',dtFileName)
    
    if (~isstruct(varargin{1}) || ischar(varargin{1})) && ~iscell(varargin{1}) 
      % A dt structue was passed.
      eigenvals = fefgGet(fg,'eigenvals',varargin{1});
    else
      % We assume that eigenvals were passed already
      eigenvals = varargin{1};
    end
    
    % This was the originial formulation described in:
    % C-F. Westin, S. Peled, H. Gubbjartsson, R. Kikinis, and F.A. Jolesz.
    % Geometrical diffusion measures for MRI from tensor basis analysis.
    % In Proceedings 5th Annual ISMRM, 1997.
    nFibers = fefgGet(fg,'nFibers');
    val     = cell(1,nFibers);

    for ii = 1:nFibers
      [val{ii}.linearity, val{ii}.planarity] = dtiComputeWestinShapes(eigenvals{ii});
    end
    
  case 'ntotalnodes'
     val = size(horzcat(fg.fibers{:}),2);      
    
  otherwise
    error('Unknown fg parameter: "%s"\n',param);
end

return
