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
    if nvararg > 1
        error('Too many input arguments.');
    end

    X = zeros(m);
    idx = 1;
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx};
    end
    f = objfun(X);


    %% Iterate:
    doagain = true; it = 0;
    while doagain
	it = it + 1;
        % Gradients and Hessians:
        tmp = zeros(m);
        for i = 1:n     % sum up
            tmp = tmp + convfft2( A(:,:,i), convfft2(A(:,:,i), X) - Y(:,:,i), 1 );
        end
        gx = tmp(:) + lambda * X(:)./sqrt(mu^2 + X(:).^2);

        % Check conditions to repeat iteration:
        if ~alphatoolow
            X = X_new;
            f = f_new;
        end
        doagain = norm(Hfun(xDelta(:))) > EPSILON && ~alphatoolow && (it < MAXIT);

    end

    % Return solution:
    Xsol.X = X;
    Xsol.W = X;         % dummy variable for compatibility with pdNCG.
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
