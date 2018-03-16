function [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu, varargin )
%XSOLVE_FISTA   Solve for X using FISTA method
%   - Core usage:
%       [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Xsolve_FISTA( ... , Xinit, Xpos, getbias )
%       Xinit:      initial value for X
%       Xpos:       constrain X to be a positive solution
%       getbias:    extract constant bias as well as X
%

    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    addpath([fpath '/helpers']);
    load([fpath '/../config/Xsolve_config.mat']); %#ok<*LOAD>
    g = huber(mu);

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end
  

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(m); b = zeros(n,1);
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx}.X;
        b = varargin{idx}.b;
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end
    
    idx = 3; getbias = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        getbias = varargin{idx};
    end

    
    
    %% Iterate:    
    t=1; W = X; u = b;
    costs = NaN(MAXIT,2);
    doagain = true;  it = 0;  count = 0;
    while doagain
	it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(m); grad_fu = zeros(n,1); R_A = zeros(m);
        for i = 1:n     % sum up
            Ri = convfft2(A(:,:,i), W) + u(i) - Y(:,:,i);
            grad_fW = grad_fW + convfft2( A(:,:,i), Ri, 1 );
            grad_fu(i) = sum(Ri(:));
            R_A = R_A + abs(fft2(A(:,:,i),m(1),m(2))).^2;
        end

        % FISTA update
        L = max(R_A(:));
        %size(W);
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        %size(X_);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X);
        if getbias
            b_ = u - grad_fu/(2*prod(m)*sqrt(n));  
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        %TODO Check conditions to repeat iteration:
        f = 0;
        for i = 1:n
            f = f + norm(convfft2(A(:,:,i), reshape(X, m)) + b(i) - Y(:,:,i), 'fro')^2/2;
        end
        costs(it,1) = f;
        costs(it,2) = g.cost(X, lambda);

        tmp = grad_fW;
        for i = 1:n
            %tmp(:,:,i) = tmp(:,:,i) + grad_fu(i); 
            tmp = tmp + grad_fu(i); %testing
        end
        delta = g.diffsubg(X, -tmp, lambda, xpos);
        delta = norm(delta(:))/sqrt(prod(m));
        if delta < EPSILON
            count = count+1;
        else
            count = 0;
        end
        doagain = count < 10 && (it < MAXIT);
    end
    

    % Return solution:
    Xsol.X = X;
    Xsol.b = b;
    Xsol.W = W;         % dummy variable for compatibility with pdNCG.
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end
