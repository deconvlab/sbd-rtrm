function [ Aout, Xsol, extras ] = Asolve_Manopt( Y, Ain, lambda, mu, varargin)
%ASolve_MANOPT     BD using Manopt solvers.
%   - Core usage:
%       [ Aout, Xsol, Stats ] = Asolve_Manopt( Y, Ain, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Asolve_Manopt( ... , Xinit, dispfun )
%

    load([fileparts(mfilename('fullpath')) '\..\config\Asolve_config.mat']);
    
    k = size(Ain); 
    if (numel(k) > 2)
        n = k(3); k = k(1:2);
    else
        n = 1;
    end
    
    Ain = Ain/norm(Ain(:));
    
    %% Handle the extra variables:    
    nvarargin = numel(varargin);
    if nvarargin > 2
        error('Too many input arguments.');
    end
    
    idx = 1;
    if nvarargin < idx || isempty(varargin{idx})
        xinit = Xsolve_pdNCG(Y, Ain, lambda, mu, []);
    else
        xinit = varargin{idx};
    end
    
    idx = 2;
    if nvarargin < idx || isempty(varargin{idx})
        dispfun = @(a, X) 0;
    else
        dispfun = varargin{idx};
    end
    
    %% Set up the problem structure for Manopt and solve
    % The package containing supplement information for the cost, egrad,
    % and ehess functions:
    suppack = var2struct(Y, k, n, lambda, mu, xinit, saveiterates);
    
    problem.M = spherefactory(prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, suppack);
    problem.egrad = @(a, store) egradfun(a, store, suppack);
    problem.ehess = @(a, u, store) ehessfun(a, u, store, suppack);
    
    options.statsfun = @(problem, a, stats, store) statsfun( problem, a, stats, store, suppack, dispfun);
    %options.stopfun = @(problem, x, info, last) stopfun(problem, x, info, last, TRTOL);
    
    % Run Manopt solver:
    [Aout, extras.cost, info, extras.options] = ManoptSolver(problem, Ain(:), options);
    
    % Produce final output:
    Aout = reshape(Aout, [k n]);
    if saveiterates
        extras.Aiter = arrayfun(@(i) info(i).A, 1:numel(info), 'UniformOutput', false);
        niter = numel(extras.Aiter);
        extras.Aiter = cell2mat(reshape(extras.Aiter, [1 1 niter]));
        extras.Xiter = arrayfun(@(i) info(i).X, 1:numel(info), 'UniformOutput', false);
        extras.Xiter = cell2mat(reshape(extras.Xiter, [1 1 niter]));
        
        Xsol.X = extras.Xiter(:,:,end);
        Xsol.W = info(end).W;
    else
        Xsol = Xsolve_pdNCG( Y, Aout, lambda, mu, suppack.xinit );
    end
end

function [ cost, store ] = costfun( a, store, suppack )
    if ~isfield(store, 'X')
        store = computeX( a, store, suppack );
    end
    
    cost = store.cost;
end

function [ egrad, store ] = egradfun( a, store, suppack )
    [Y, k, n] = unpack( suppack, 'Y', 'k', 'n' );
    if ~isfield(store, 'X')
        store = computeX( a, store, suppack );
    end
    
    m = size(store.X);
    egrad = zeros(prod(k)*n,1);
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        tmp = convfft2( store.X, convfft2( reshape(a(idx), k), store.X ) - Y(:,:,i), 1, m+k-1, m);
        tmp = tmp(1:k(1), 1:k(2));
        egrad(idx) = tmp(:);
    end
end

function [ ehess, store ] = ehessfun( a, u, store, suppack )
    [Y, k, n, lambda, mu] = unpack( suppack, 'Y', 'k', 'n', 'lambda', 'mu');
    if ~isfield(store, 'X')
        store = computeX( a, store, suppack );
    end

    ehess = H_function( u, Y, reshape(a, [k n]), store.X, lambda, mu );
end

function [ store ] = computeX( a, store, suppack )
    % Updates the cache to store X*(A), and the active-set whenever a new
    % a new iteration by the trust-region method needs it.
    
    [Y, k, n, lambda, mu, xinit] = unpack( suppack, ...
        'Y', 'k', 'n', 'lambda', 'mu', 'xinit');
    sol = Xsolve_pdNCG( Y, reshape(a, [k n]), lambda, mu, xinit );
    
    store.X = sol.X;
    store.W = sol.W;
    store.cost = sol.f;
end

function [ stats ] = statsfun( problem, a, stats, store, suppack, dispfun) %#ok<INUSL>
    if suppack.saveiterates
        stats.A = reshape(a, [suppack.k suppack.n]);
        stats.X = store.X;      % So X could be returned at the end.
        stats.W = store.W;
    end
    dispfun(a, store.X);
end

function s = var2struct(varargin)
  names = arrayfun(@inputname,1:nargin,'UniformOutput',false);
  s = cell2struct(varargin,names,2);
end
