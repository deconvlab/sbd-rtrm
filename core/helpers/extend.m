function [ out ] = extend( in, outsz, varargin )
%EXTEND     The extension map and its adjoint.
    if numel(varargin) > 1
        error('Too many input arguments.');
    end
    if numel(varargin) == 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    
    if adj == false
        if sum(outsz < size(in))
            error('OUTSZ needs to be geq than SIZE(IN) if ''ADJ == FALSE''.');
        end
        out = zeros(outsz);
        out(1:size(in,1), 1:size(in,2)) = in;
    else
        if sum(outsz > size(in))
            error('OUTSZ needs to be leq than SIZE(IN) if ''ADJ == TRUE''.');
        end
        out = in(1:outsz(1), 1:outsz(2));
    end

end

