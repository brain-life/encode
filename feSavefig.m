function figDir = feSavefig(h,varargin)
% Saves a figure for publishing purpose.
%
%  function figDir = savefig(h,varargin)
% 
% INPUTS: Must be pairs, e.g., 'name',value
%  
%   figName  - the name of the figure file. 
%              Default: sprintf('feFig_%s',get(h,'Name'))
%   figDir   - Name of the subfolder where to save this figure.
%              Default: current dir ('.')
%   figType  - fig, eps, jpg, png
%              Default: 'png' fastest smallest file, low resolution.
%   verbose  - display on screen the print command being invoked.
% 
% NOTES:
%   This function invokes a series of print commands:
%   print(h, '-cmyk', '-painters','-depsc2','-tiff','-r500', '-noui', 'figureName')
%
% EXAMPLE:
%   feSavefig(figureHandle,'verbose','yes','figName','myfig','figDir','/path/to/fig/folder/');
%
%
% Copyright (2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% set up default variables:
figName           = sprintf('feFig_%s',get(h,'Name')); % the name of the figure file
figDir            = '.'; % the subfolder where to save this figure
figType           = 'png';
verbose           = 'yes'; % 'yes', 'no', display on screen what it is going on

if ~isempty(varargin)
  if mod(length(varargin),2), error('varargin must be pairs'); end
  for ii=1:2:(length(varargin)-1)
    eval(sprintf('%s = ''%s'';',varargin{ii}, varargin{ii+1}));
  end
end

% make the figure dir if it does not exist:
if ~isdir(figDir), mkdir(figDir);end

% Create a print command that will save the figure
switch figType
  case {'png'}
    printCommand = ...
    sprintf('print(%s, ''-painters'',''-dpng'', ''-noui'', ''%s'')', ...
    num2str(h),fullfile(figDir,figName));
  
  case {'jpg'}
    printCommand = ...
    sprintf('print(%s, ''-djpeg95'',''-r500'', ''-noui'', ''%s'')', ...
    num2str(h),fullfile(figDir,figName));

  case {'eps'}
    printCommand = ...
    sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', ...
    num2str(h),fullfile(figDir,figName));
  
  case {'fig'}
    printCommand = ...
    sprintf('saveas(%s, ''%s'', ''fig'')', num2str(h),figName);
    
  otherwise
    error('[%s] Cannot save figure with type set to: %s', mfilename,figType)
end

if strcmpi(verbose,'yes')
 fprintf('[%s] Saving eps figure %s/%s\nUsing command: \n%s...\n', ...
 mfilename,figDir,figName,printCommand);
end

% Do the printing here:
eval(printCommand);

% Delete output if it was nto requested
if (nargout < 1), clear figDir;end
return
