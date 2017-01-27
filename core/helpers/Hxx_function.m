function [ out ] = Hxx_function( v, m, A, Hdiag )
    %% 1. Apply rho'*Ca'*P'*P*Ca*rho to v
    % Reshape and extend v
    k = [size(A,1) size(A,2)]; n = size(A,3);
    tmpsz = m + k - 1;
    rhov = extend(reshape(v, m), tmpsz);
    
    % For each slice, apply the rest of the operations
    tmp = zeros(m);
    delta = tmpsz - m;
    mask = false(tmpsz); 
    mask(floor(delta(1)/2)+(1:m(1)), floor(delta(2)/2)+(1:m(2))) = true;
    for i = 1:n
        % apply P'*P*Ca to rhov
        currslice = mask.*cconvfft2(A(:,:,i),rhov);
        
        % apply rho'*Ca' and save
        currslice = extend(cconvfft2(A(:,:,i),currslice,tmpsz,'left'),m,1);
        tmp = tmp + currslice;
    end
    out = tmp(:) + Hdiag.*v;
end