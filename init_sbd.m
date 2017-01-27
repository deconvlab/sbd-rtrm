function [] = init_sbd( mode, setdefconfig )
%INIT_SBD   Initializes subdirectories and default config settings.
%
%   Optional arguments:
%       mode - to turn off messages, set to 'quiet'. 
%           Verbose by default.
%
%       setdefconfig - whether to apply default config settings. 
%           Default: true.
%

    if nargin < 1;  mode = 'verbose';       end;
    if nargin < 2;  setdefconfig = true;    end;
    
    % Check if ManOpt has been imported
    if exist('trustregions', 'file') ~= 2
        error('Error initializing SBD package.\nManopt hasn''t been imported yet!  http://www.manopt.org%s','');
    end
    
    % Add subdirectories to path
    fp = [fileparts(mfilename('fullpath')) '\'];
    addpath(fp);
    for d = {'core', 'utils', 'config'}
        addpath(genpath([fp d{1}]));
    end
    
    % Apply default config settings
    if setdefconfig;    default_config_settings(mode);  end;
    
    if ~strcmp(mode, 'quiet')
        disp('Subdirectories and config settings initialized.');
    end
end