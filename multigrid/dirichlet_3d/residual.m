function r = residual(u,f,h)
%RESIDUAL
% This function computes the difference between the right-hand side f and
% the Laplacian of u.
%
% Syntax: r = residual(u,f,h)
%
% Inputs:
%   u: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   h: The grid spacing. Scalar.
%
% Outputs:
%   r: The residual 
%
% See also: LAPLACIAN
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 31st of October 2022
    assert(sum(size(u) == size(f)) == 3,...
        "The shape of the iterate u and the right-hand side f must match");
    r = f-laplacian(u,f,h);
end