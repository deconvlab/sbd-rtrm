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
    g = huber(mu);

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 2
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(m);
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx};
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end


    %% Iterate:
    t=1; W = X;
    costs = NaN(MAXIT,1);
    doagain = true;  it = 0;  count = 0;
    while doagain
	it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(m); R_A = zeros(m);
        for i = 1:n     % sum up
            grad_fW = grad_fW + convfft2( A(:,:,i), convfft2(A(:,:,i), W) - Y(:,:,i), 1 );
            R_A = R_A + fft2(A(:,:,i),m(1),m(2));
        end

        % FISTA update
        L = max(R_A(:));
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X)
        X = X_; t = t_;

        %TODO Check conditions to repeat iteration:
        f = 0;
        for i = 1:n
            f = f + norm(convfft2(A(:,:,i), reshape(X, m)) - Y(:,:,i), 'fro')^2/2;
        end
        f = f + g.cost(X, lambda);
        costs(it) = f;

        delta = g.diffsubg(X, grad_fW, lambda, xpos);
        if norm(delta(:))/sqrt(prox(m)) < EPSILON
            count = count+1;
        else
            count = 0;
        end
        doagain = count > 10 && (it < MAXIT);
    end

    % Return solution:
    Xsol.X = X;
    Xsol.W = W;         % dummy variable for compatibility with pdNCG.
    Xsol.f = f;
    info.numit = it;
    info.costs = costs(1:it);
end
