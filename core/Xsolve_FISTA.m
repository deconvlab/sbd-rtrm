function [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu, varargin )
%XSOLVE_FISTA   Solve for X using FISTA method
%   - Core usage:
%       [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Xsolve_FISTA( ... , Xinit, Xpos )
%       Xinit: initial value for X
%       Xpos: constrain X to be a positive solution
%

    % Initialize variables and function handles:
    addpath('helpers')
    load([fileparts(mfilename('fullpath')) '\..\config\Xsolve_config.mat']);

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    objfun = @(X) obj_function ( X, A, Y, lambda, mu );

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 2
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(m);
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx};
    end
    f = objfun(X);

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end


    %% Iterate:
    doagain = true; it = 0;
    t=1; W = X;
    while doagain
	it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(m);
        for i = 1:n     % sum up
            grad_fW = grad_fW + convfft2( A(:,:,i), convfft2(A(:,:,i), W) - Y(:,:,i), 1 );
        end

        % FISTA update
        L = 1;  %TODO Work out the Lipschitz constant!
        X_ = prox_huber(W - 1/L*grad_fW, lambda/L, mu, xpos);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X)
        X = X_; t = t_;

        %TODO Check conditions to repeat iteration:
        if ~alphatoolow
            X = X_new;
            f = f_new;
        end
        doagain = norm(Hfun(xDelta(:))) > EPSILON && ~alphatoolow && (it < MAXIT);
    end

    % Return solution:
    Xsol.X = X;
    Xsol.W = W;         % dummy variable for compatibility with pdNCG.
    Xsol.f = f;
    info.numit = it;
end

function [ out ] = obj_function ( X, A, Y, lambda, mu )
    m = size(Y);

    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    out = 0;
    for i = 1:n
        out = out + norm(convfft2(A(:,:,i), reshape(X, m)) - Y(:,:,i), 'fro')^2/2;
    end
    out = out + lambda*huber(X);
end
