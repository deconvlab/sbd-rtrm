function [ h ] = huber( mu )
    if nargin < 1 || isempty(mu)
        mu = 1e-6;
    end

    h.cost = @(x, lambda) cost(x, lambda, mu);
    h.prox = @(x, lambda, xpos) prox (x, lambda, mu, xpos);
    h.diffsubg = @(x, y, lambda, xpos) diffsubg(x, y, lambda, mu, xpos);
end

function [ hx ] = cost(x, lambda, mu)
    leq = abs(x) <= mu;
    hx = sum(x(leq).^2/2/mu) + sum(abs(x(~leq)) - mu/2);
    hx = lambda*hx;
end

function [ proxh ] = prox(x, lambda, mu, xpos)
    if nargin < 4 || isempty(xpos)
        xpos = false;
    end

    leq = abs(x) <= lambda + mu;
    proxh = zeros(size(x));
    proxh(leq) = x(leq)./(1+lambda/mu);
    proxh(~leq) = x(~leq) - lambda*sign(x(~leq));

    if xpos;  proxh = max(proxh,0);  end
end

function [ eps ] = diffsubg(x, y, lambda, mu, xpos)
%DIFFSUBG  Difference from y to the subgradient of the loss function at x.
%  In this case the Huber loss function is smooth so only need to subtract
%  from the gradient.

    leq = abs(x) <= mu;
    subg = zeros(size(x));

    subg(leq) = x(leq)/mu;
    subg(~leq) = sign(x(~leq));
    if xpos;  subg(x<0) = 0;  end

    eps = y - lambda*subg;
end