function [ proxh ] = prox_huber(x, lambda, mu, xpos)
    if nargin < 3 || isempty(mu)
        mu = 1e-6;
    end
    if nargin < 4 || isempty(xpos)
        xpos = false;
    end

    leq = abs(x) <= lambda + mu;
    proxh = zeros(size(x));
    proxh(leq) = x./(1+lambda/mu);
    proxh(~leq) = x - lambda*sign(x);

    if xpos;  proxh = max(proxh,0);  end
end