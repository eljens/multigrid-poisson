function uc = restrict(uf)
%RESTRICT
% This function selects every second point in a domain in each dimension.
%
% Syntax: uc = restrict(uf)
%
% Inputs:
%   uf: The current iterate of the solution on the finer domain. 3D matrix.
%
% Outputs:
%   uc: The restriction of uf on the coarse domain. 3D Matrix.
%
% See also: INTERPOLATE
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    if (~isscalar(uf))
        assert(mod(size(uf,1),2) == 1,"First dimension of of must be 2^l+1")
        assert(mod(size(uf,2),2) == 1,"Second dimension of of must be 2^l+1")
        assert(mod(size(uf,3),2) == 1,"Third dimension of of must be 2^l+1")
    end
    uc = uf(1:2:end,1:2:end,1:2:end);
end