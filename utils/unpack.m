function [ varargout ] = unpack( s, varargin )
%UNPACK     Unpack structure fields into variables.
%   - Usage:
%       [ out1, out2, ... ] = unpack( s, fieldname1, fieldname2, ... )
%   
%   - e.g.: 
%       s.A = 1;    s.B = 'foo';    s.C = true;
%       [out1, out2, out3] = unpack( s, 'B', 'A', 'D' );
%       
%     Then,
%       out1 = 'foo',   out2 = 1,   out3 = [].
%

    n = numel(varargin);
    varargout = cell(n,1);
    
    for i = 1:numel(varargin)
        currvar = varargin{i};
        if ~ischar(currvar)
            error('Field names should all be strings');
        elseif isfield(s, currvar)
            varargout{i} = s.(currvar);
        else
            varargout{i} = [];
        end
    end

end

