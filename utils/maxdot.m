function [ max_dot ] = maxdot( A0, Ahat, varargin )
%MAXDOT   Max. abs. dot prod between A0 and Ahat over shifts.
%
%   Usage:      [ max_dot ] = maxdot(A0, Ahat)
%   
%   Optional input:      [...] = maxdot(A0, Ahat, kplus)
%       Allows user to adjust the amount of allowable shifts.
%
    k = size(A0); 
    if (numel(k) > 2)
        n = k(3); k = k(1:2);
    else
        n = 1;
    end
    
    nvarargin = numel(varargin);
    kplus = ceil(0.5*k);
    if (nvarargin >= 1)
        kplus = varargin{1};
    end
    
    Apad = zeros([k+2*kplus n]);
    Apad(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :) = Ahat;
    
    absdotprod_mtx = zeros(2*kplus + 1);
    for i = -kplus(1):kplus(1)
        for j = -kplus(2):kplus(2)
            idx = [i j]+kplus+1;
            A_try = Apad(idx(1)+(0:k(1)-1),idx(2)+(0:k(2)-1),:);

            tempvar = abs(dot( A0(:), A_try(:) ))/(norm(A0(:))*norm(A_try(:)));
            absdotprod_mtx(idx(1),idx(2)) = tempvar;
        end
    end
    max_dot = max(absdotprod_mtx(:));
end

