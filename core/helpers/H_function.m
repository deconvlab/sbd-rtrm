function [ H_v ] = H_function( v, Y, A, X, lambda, mu )
%H_FUNCTION     Apply the Euclidean Hessian to a vector.
    load([fileparts(mfilename('fullpath')) '/../../config/Hfunction_config.mat']);
    
    m = size(Y); k = size(A); n = size(Y,3);
    m = m(1:2); k = k(1:2);
    tmpsz = m + k - 1;
    
    Haa_v = zeros([k n]);
    r = zeros([m n]);
    Hxa_v = zeros(m);
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        vi = reshape(v(idx), k);
        Haa_v(:,:,i) = convfft2( X, convfft2(X, vi, 0, tmpsz, m), 1, tmpsz, k);
    
        r(:,:,i) = convfft2(A(:,:,i),X) - Y(:,:,i);
        Hxa_v = Hxa_v ...
            + convfft2( A(:,:,i), convfft2(X, vi, 0, tmpsz, m), 1) ...
            + Hxres( r(:,:,i), vi, k, m );
    end
    
    hesspendiag = lambda * mu^2*(mu^2 + X(:).^2).^(-3/2);    
    Hxxfun = @(u) Hxx_function(u, m, A, hesspendiag);
    pcgprecond = @(u) u./(1 + hesspendiag);
    [HxxInv_Hxa_v,~] = pcg(Hxxfun, Hxa_v(:), PCGTOL, PCGIT, pcgprecond);
    HxxInv_Hxa_v = reshape(HxxInv_Hxa_v, m);
    
    % Hax = Ikm'*(CA*CX' + (CA*CX-CY)*Pi);
    
    Hax_HxxInv_Hxa_v = zeros([k n]);
    for i = 1:n
        Hax_HxxInv_Hxa_v(:,:,i) = ...
            convfft2( X, convfft2(A(:,:,i), HxxInv_Hxa_v, 0, tmpsz, m), 1, tmpsz, k) ...
            + Hxres( r(:,:,i), HxxInv_Hxa_v, m, k );
    end
    
    H_v = Haa_v - Hax_HxxInv_Hxa_v;
    H_v = H_v(:);

end

function [ out ] = Hxres( res, in, insz, outsz )
    tmpsz = outsz + insz - 1;
    Vhat = fft2(extend(reshape(in, insz), tmpsz));
    reshat = fft2(shrink(res, tmpsz, 1));
    out = extend(ifft2(reshat.*conj(Vhat)), outsz, 1 );
end