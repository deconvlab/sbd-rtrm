function [ Aout, Xsol, extras ] = Asolve_Manopt( Y, Ain, lambda, varargin )
%ASolve_MANOPT     BD using Manopt solvers.
%   - Core usage:
%       [ Aout, Xsol, Stats ] = Asolve_Manopt( Y, Ain, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Asolve_Manopt( ... , Xinit, Xpos, dispfun )
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
        xpos = false;
    else
        xpos = varargin{idx};
    end

    idx = 3;
    if nvarargin < idx || isempty(varargin{idx})
        dispfun = @(a, X) 0;
    else
        dispfun = varargin{idx};
    end

    %% Set up the problem structure for Manopt and solve
    % The package containing supplement information for the cost, egrad,
    % and ehess functions:
    %{
    suppack.Y = Y;
    suppack.k = k;
    suppack.n = n;
    suppack.mu = mu;
    suppack.xinit = xinit;
    suppack.saveiterates = saveiterates;
    %}

    problem.M = spherefactory(prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, Y, k, n, lambda, mu, xinit, xpos );
    problem.egrad = @(a, store) egradfun(a, store,  Y, k, n, lambda, mu, xinit, xpos );
    problem.ehess = @(a, u, store) ehessfun(a, u, store,  Y, k, n, lambda, mu, xinit, xpos );

    options.statsfun = @(problem, a, stats, store) statsfun( problem, a, stats, store, k, n, saveiterates, dispfun);
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
        Xsol = Xsolve_FISTA( Y, Aout, lambda, mu, suppack.xinit, xpos );
    end
end

function [ cost, store ] = costfun( a, store, Y, k, n, lambda, mu, xinit, xpos )
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos );
    end

    cost = store.cost;
end

function [ egrad, store ] = egradfun( a, store, Y, k, n, lambda, mu, xinit, xpos )
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos );
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

function [ ehess, store ] = ehessfun( a, u, store, Y, k, n, lambda, mu, xinit, xpos )
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos );
    end

    ehess = H_function( u, Y, reshape(a, [k n]), store.X, lambda, mu );
end

function [ store ] = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos )
    % Updates the cache to store X*(A), and the active-set whenever a new
    % a new iteration by the trust-region method needs it.

    sol = Xsolve_pdNCG( Y, reshape(a, [k n]), lambda, mu, xinit, xpos );

    store.X = sol.X;
    store.W = sol.W;
    store.cost = sol.f;
end

function [ stats ] = statsfun( problem, a, stats, store, k, n, saveiterates, dispfun) %#ok<INUSL>
    if saveiterates
        stats.A = reshape(a, [k n]);
        stats.X = store.X;      % So X could be returned at the end.
        stats.W = store.W;
    end
    dispfun(a, store.X);
end