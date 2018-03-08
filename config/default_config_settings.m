function [] = default_config_settings( mode )   
%DEFAULT_CONFIG_SETTINGS    Default settings used in core BD2 functions.
%   Change the default configuration settings for various solvers here.
%
%   Optional arguments:
%       mode - to turn off messages, set to 'quiet'. 
%           Verbose by default.
%

    fp = [fileparts(mfilename('fullpath')) '/']; %/ for linux, \ for windows.
    clean = @() eval('clearvars -except fp mode clean shared');
    
    %% Shared settings
    shared.mu = 1e-6;           % pseudo-Huber approximation parameter
    
    %% Settings for ASOLVE_MANOPT.M
    
    mu = shared.mu;	
    saveiterates = false;       % keep Asolve iterates?
    
    % Manopt settings:
    options.verbosity = 2;
    options.tolgradnorm = 1e-8;
    options.maxiter = 3;          % maximum iterations for RTRM
    options.storedepth = 5;

    %   for gradient methods:
    options.linesearch = @linesearch;
    options.ls_contraction_factor = 0.2;
    options.ls_suff_decr = 1e-3; %#ok<*STRNU>
    
    %   descent method:
    ManoptSolver = @trustregions;
    
    save([fp 'Asolve_config.mat']); clean;

    %% Settings for XSOLVE_PDNCG.M
    
    EPSILON = 1e-8;    % Tolerance to stop the x-solver.
    ALPHATOL = 1e-10;  % When alpha gets too small, stop.
    C2 = 1e-6;         % How much the obj. should decrease; 0 < C2 < 0.5.
    C3 = 2e-1;         % Rate of decrease in alpha.
    PCGTOL = 1e-8;
    PCGIT = 2e2;       % Maximum iterations for PCG
    MAXIT = 1e2;       % Maximum iterations for XSolve

    save([fp 'Xsolve_config.mat']); clean;
    
    %% Settings for H_FUNCTION.M
    
    PCGTOL = 1e-8;
    PCGIT = 2e2;   % Maximum iterations for PCG

    save([fp 'Hfunction_config.mat']); clean;
    
    %% Confirmation:
    if ~strcmp(mode, 'quiet')
        disp('Default configuration settings applied.')
    end
    
end

%#ok<*NASGU>
%#ok<*VUNUS>
