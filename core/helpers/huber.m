function [ hx ] = huber( x, mu )
    if nargin < 2 || isempty(mu)
        mu = 1e-6;
    end
    leq = abs(x) <= mu;
    hx = sum(x(leq).^2/2/mu) + sum(abs(x(~leq)) - mu/2);
end