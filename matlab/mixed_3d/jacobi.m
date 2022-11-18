function unew = jacobi(u,f,h,gx1,gxn,gy1,gyn)
%JACOBI
% This function performs one Jacobi iteration for a uniform discretization
% of the 3D Poisson problem on the form $\Delta u = f$. The function uses
% weight omega=2/3.
%
% Syntax: unew = jacobi(u,f,h)
%
% Inputs:
%   u:      The current iterate of the solution. 3D matrix.
%   f:      The right hand side. 3D matrix.
%   h:      The grid spacing. Scalar.
%   gx1:    The Neumann condition on the west boundary. 2D matrix.
%   gxn:    The Neumann condition on the east boundary. 2D matrix.
%   gy1:    The Neumann condition on the south boundary. 2D matrix.
%   gyn:    The Neumann condition on the north boundary. 2D matrix.
%
% Outputs:
%   unew:   The next iterate of the solution.
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

    % Western Neumann condition
    unew(1,j,k) = (1-omega)*u(1,j,k) + omega*(1/6)*(2*u(2,j,k)...
        +u(1,j-1,k)+u(1,j+1,k)+u(1,j,k-1)+u(1,j,k+1)...
        -f(1,j,k)*hsq-2*h*gx1(1,j,k));

    % Eastern Neumann condition
    unew(end,j,k) = (1-omega)*u(end,j,k) + omega*(1/6)*(2*u(end-1,j,k)...
      +u(end,j-1,k)+u(end,j+1,k)+u(end,j,k-1)+u(end,j,k+1)...
      -f(end,j,k)*hsq+2*h*gxn(1,j,k));

    % Southern Neumann condition
    unew(i,1,k) = (1-omega)*u(i,1,k) + omega*(1/6)*(u(i-1,1,k)+u(i+1,1,k)...
        +2*u(i,2,k)+u(i,1,k-1)+u(i,1,k+1)...
        -f(i,1,k)*hsq-2*h*gy1(i,1,k));

    % Northern Neumann condition
    unew(i,end,k) = (1-omega)*u(i,end,k) + omega*(1/6)*(u(i-1,end,k)+u(i+1,end,k)...
        +2*u(i,end-1,k)+u(i,end,k-1)+u(i,end,k+1)...
        -f(i,end,k)*hsq+2*h*gyn(i,1,k));

    % South west Neumann condition
    unew(1,1,k) = (1-omega)*u(1,1,k) + omega*(1/6)*(2*u(2,1,k)...
        +2*u(1,2,k)+u(1,1,k-1)+u(1,1,k+1)...
        -f(1,1,k)*hsq-2*h*gx1(1,1,k)-2*h*gy1(1,1,k));

    % South east Neumann condition
    unew(end,1,k) = (1-omega)*u(end,1,k) + omega*(1/6)*(2*u(end-1,1,k)...
        +2*u(end,2,k)+u(end,1,k-1)+u(end,1,k+1)...
        -f(end,1,k)*hsq+2*h*gxn(1,1,k)-2*h*gy1(end,1,k));

    % North west Neumann condition
    unew(1,end,k) = (1-omega)*u(1,end,k) + omega*(1/6)*(2*u(2,end,k)...
        +2*u(1,end-1,k)+u(1,end,k-1)+u(1,end,k+1)...
        -f(1,end,k)*hsq-2*h*gx1(1,end,k)+2*h*gyn(1,1,k));

    % North east Neumann condition
    unew(end,end,k) = (1-omega)*u(end,end,k) + omega*(1/6)*(2*u(end-1,end,k)...
        +2*u(end,end-1,k)+u(end,end,k-1)+u(end,end,k+1)...
        -f(end,end,k)*hsq+2*h*gxn(1,end,k)+2*h*gyn(end,1,k));
end