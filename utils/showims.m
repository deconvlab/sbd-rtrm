function showims( Y, A0, X0, A, X, k, kplus, idx )
%SHOWIMS    Show images after each iteration.
    A = reshape(A, [k size(Y,3)]);

    Y_hat = convfft2(A(:,:,idx),X);
    subplot(321); imagesc(abs(Y(:,:,idx)));     ylabel('abs(Y)'); title('Original');
    subplot(322); imagesc(abs(Y_hat));          title('Recovered');
    
    subplot(324); imagesc(abs(A(:,:,idx)));
    subplot(323); imagesc(abs(A0(:,:,idx)));
    
    if ~isempty(kplus)      % i.e. we're in Phase II
        ylabel('abs(A)  (lifted)'); 
        X_hat = circshift(X,kplus);
    else
        ylabel('abs(A)'); 
        X_hat = X;
    end
    
    subplot(325); imagesc(abs(X0)); ylabel('abs(X)');
    subplot(326); imagesc(abs(X_hat));
    
    drawnow;
end

