function [ Y ] = convfft2( A, B, varargin )
%CCONVFFT2	FFT implementation of 2D linear convolution
%   Only non-extended portion of convolution (size of B) kept.
%
%   Y = convfft2(A, B)	convolves A with B.
%   In 1D vector notation,      c = P*C{iota*a}*(rho*b), 
%   where rho extends B to size(A)+size(B)-1, and P shrinks the 
%   window observed from the output of Ca*rho*b to SIZE(B).
%
%   Y = convfft2(..., adj, tmpsz, outsz)  
%   applies the adjoint operation if ADJ==TRUE: rho'*C{iota*a}'*P'.
%   TMPSZ and OUTSZ sets the extension and output sizes manually.
%

    if numel(varargin) > 3
        error('Too many input arguments.');
    end
    if numel(varargin) >= 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        tmpsz = varargin{2};
    else
        tmpsz = size(B) + size(A) - 1;
    end
    if numel(varargin) >= 3 && ~isempty(varargin{3})
        outsz = varargin{3};
    else
        outsz = size(B);
    end
    
    if adj == false
        Y = shrink(cconvfft2(A, extend(B,tmpsz)), outsz);
    else
        Y = extend(cconvfft2(A, shrink(B, tmpsz, 1), tmpsz, 'left'), outsz, 1);
    end
end





