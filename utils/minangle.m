function [ minang ] = minangle( A0, Ahat, varargin )
%MINANGLE   Min. angle between A0 and Ahat over shifts.
%
%   Usage:      [ minang ] = minangle(A0, Ahat)
%   
%   Optional input:      [...] = minangle(A0, Ahat, kplus)
%       Allows user to adjust the amount of allowable shifts.
%
    
    if isempty(A0)
        minang = NaN;
        return;
    end

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
    
    mang_mtx = zeros(2*kplus + 1);
    for i = -kplus(1):kplus(1)
        for j = -kplus(2):kplus(2)
            idx = [i j]+kplus+1;
            A_try = Apad(idx(1)+(0:k(1)-1),idx(2)+(0:k(2)-1),:);

            tempvar = abs(dot( A0(:), A_try(:) ))/(norm(A0(:))*norm(A_try(:)));
            mang_mtx(idx(1),idx(2)) = acos(tempvar);
        end
    end
    minang = min(mang_mtx(:));
end

