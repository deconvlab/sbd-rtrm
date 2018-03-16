function [ Xsol, info ] = Xsolve_pdNCG( Y, A, lambda, mu, varargin )
%XSOLVE2_PDNCG	Solves for X*(A) by using a primal-dual Newton CG method.
%   - Core usage:
%       [ Xsol, info ] = Xsolve_pdNCG( Y, A, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Xsolve_pdNCG( ... , Xinit, xpos, getbias)
%
%   Algorithm from (Fountoulakis and Gondzio '14).

    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    addpath([fpath '/helpers']);
    load([fpath '/../config/Xsolve_config.mat']); %#ok<*LOAD>

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    objfun = @(X) obj_function ( X, A, Y, lambda, mu );

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    X = zeros(m); W = zeros(m);
    idx = 1;
    if nvararg >= idx && ~isempty(varargin{idx})
        if isfield(varargin{idx}, 'X') && ~isempty(varargin{idx}.X)
            X = varargin{idx}.X;
        end
        if isfield(varargin{idx}, 'W') && ~isempty(varargin{idx}.W)
            W = varargin{idx}.W;
        end
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

        D = 1./sqrt(mu^2 + X(:).^2);
        Hdiag = lambda*D.*(1 - D.*X(:).*W(:));
        Hfun = @(v) Hxx_function(v, m, A, Hdiag);
        PCGPRECOND = @(v) v./(Hdiag + 1);

        % Solve for xDelta using PCG:
        [xDelta,~] = pcg(Hfun, -gx, PCGTOL, PCGIT, PCGPRECOND);
        xDelta = reshape(xDelta, m);

        % Update the dual variable:
        wDelta = D.*( 1 - D.*X(:).*W(:) ).*xDelta(:) - ( W(:) - D.*X(:) );
        W = W + reshape(wDelta, m);
        W = min(abs(W), 1).*sign(W);

        % Update the primal variable by backtracking:
        alpha = 1/C3; f_new = Inf; alphatoolow = false;
        while f_new > f - C2*alpha*norm(Hfun(xDelta(:)))^2 && ~alphatoolow
            alpha = C3*alpha;
            X_new = X + alpha*xDelta;
            f_new = objfun(X_new);

            % [f_new f C2*alpha*norm(Hfun(xDelta(:)))^2]
            alphatoolow = alpha < ALPHATOL;
        end

        % Check conditions to repeat iteration:
        if ~alphatoolow
            X = X_new;
            f = f_new;
        end
        doagain = norm(Hfun(xDelta(:))) > EPSILON && ~alphatoolow && (it < MAXIT);

    end

    % Return solution:
    Xsol.X = X;
    Xsol.W = W;
    Xsol.f = f;
    Xsol.b = 0; %compatible with FISTA
    info.numit = it;
    info.alphatoolow = alphatoolow;
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
    out = out + lambda.*sum(sqrt(mu^2 + X(:).^2) - mu);
end
