function [ out ] = shrink( in, outsz, varargin )
%RESTRICT   Restriction of observation window and its adjoint.
    if numel(varargin) > 1
        error('Too many input arguments.');
    end
    if numel(varargin) == 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    
    if adj == false
        if sum(outsz > size(in))
            error('OUTSZ needs to be leq than SIZE(IN) if ''ADJ == FALSE''.');
        end
        delta = size(in) - outsz;
        out = in(floor(delta(1)/2)+(1:outsz(1)), floor(delta(2)/2)+(1:outsz(2)));
    else
        if sum(outsz < size(in))
            error('OUTSZ needs to be geq than SIZE(IN) if ''ADJ == TRUE''.');
        end
        delta = outsz - size(in);
        out = zeros(outsz);
        out(floor(delta(1)/2)+(1:size(in,1)), floor(delta(2)/2)+(1:1:size(in,2))) = in;
    end
end

