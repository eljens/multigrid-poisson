function A = laplacian(u,f,h)
%LAPLACIAN
% This function approximates the Laplacian operator on a 3D domain. This
% function imposes Diriclet boundary conditions on all surface points.
%
% Syntax: A = laplacian(u,f,h)
%
% Inputs:
%   u: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   h: The grid spacing. Scalar.
%
% Outputs:
%   A: An approximation of the Laplacian of the solution iterate u. The
%   residuals of Au=f may be computed ad r = f-Au.
%
% See also: RESIDUAL
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    assert(sum(size(u) == size(f)) == 3,...
        "The shape of the iterate u and the right-hand side f must match");
    i = 2:size(u,1)-1;
    j = 2:size(u,2)-1;
    k = 2:size(u,3)-1;
    A = f;
    A(i,j,k) = (u(i-1,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i,j+1,k)...
        +u(i,j,k-1)+u(i,j,k+1)-6*u(i,j,k))./(h^2);
end