function showimsrec( Y, A, X, k, kplus, idx )
%SHOWIMSREC     Show one recovered slice after each iteration.
    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    A = reshape(A, [k n]);
    X = reshape(X, m);

    Y_hat = convfft2(A(:,:,idx),X);
    subplot(221); imagesc(abs(Y(:,:,idx)));	title('abs(Y) - Observed');
    subplot(223); imagesc(abs(Y_hat));      title('abs(Y) - Recovered');
    
    subplot(222); imagesc(abs(A(:,:,idx)));
    
    if ~isempty(kplus)      % i.e. we're in Phase II
        title('abs(A) - Lifted'); 
        X_hat = circshift(X,kplus);
    else
        title('abs(A)'); 
        X_hat = X;
    end
    
    subplot(224); imagesc(abs(X_hat));  title('abs(X)');
    
    drawnow;
end



