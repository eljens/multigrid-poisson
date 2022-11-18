function unew = jacobi(u,f,h)
%JACOBI
% This function performs one Jacobi iteration for a uniform discretization
% of the 3D Poisson problem on the form $\Delta u = f$. The function uses
% weight $\omega=\frac{2}{3}$.
%
% Syntax: unew = jacobi(u,f,h)
%
% Inputs:
%   u: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   h: The grid spacing. Scalar.
%
% Outputs:
%   unew: The next iterate of the solution.
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    omega = 2/3;
    i = 2:size(u,1)-1;
    j = 2:size(u,2)-1;
    k = 2:size(u,3)-1;
    % The square of the grid spacing
    hsq = h^2;

    % We need a copy since the Dirichlet boundary points are not updated
    unew = u;

    % Jacobi iteration for all interior grid points
    unew(i,j,k) = (1-omega)*u(i,j,k) + omega*(1/6)*(u(i-1,j,k)+u(i+1,j,k)...
        +u(i,j-1,k)+u(i,j+1,k)+u(i,j,k-1)+u(i,j,k+1)-f(i,j,k)*hsq);
end