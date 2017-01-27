function show_realdata_results( Y, Aout, Xout, absflag )
%SHOW_REALDATA_RESULTS  Shows some results from one output slice.
%   absflag: show the absolute value of the variables - 0 or 1.
    
    m = size(Y);
    k = size(Aout);
    
    if absflag
        vfun = @(var) abs(var);
        tfun = @(str) ['abs(' str ')'];
    else
        vfun = @(var) var;
        tfun = @(str) str;
    end
    
    subplot(321);   
    imagesc(vfun(Y));  
    title(tfun('Y'));
    
    subplot(323);
    mpad = m + k - 1;
    Xpad = zeros([m 1]);
    Xpad(1:m(1), 1:m(2)) = Xout;
    Ypad = cconvfft2(Aout, Xpad);
    
    imagesc(vfun( Ypad(1:m(1), 1:m(2)) ));
    title(tfun(' conv2(A, X) '));
    
    
    subplot(324);
    imagesc(vfun(Aout));
    title(['Recovered ' tfun('A')]);

    
    subplot(322);
    imagesc(vfun(Xout));
    title(['Recovered ' tfun('X')]);
    
    
    subplot(325);
    imagesc(abs(fftshift(fft2( Y ))));
    title('abs(fft2( Y ))');
    
    
    subplot(326);
    imagesc(abs(fftshift(fft2( Aout, m(1), m(2) ))));
    title('abs(fft2( A, size(Y) ))');
end

