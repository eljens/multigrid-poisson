function A = laplacian(u,f,h,gx1,gxn,gy1,gyn)
%LAPLACIAN
% This function approximates the Laplacian operator on a 3D domain. This
% function imposes Diriclet boundary conditions on the horizontal
% boundaries (top and bottom. It imposes Neumann conditions on the vertical
% boundaries (east, west, north, and south).
%
% Syntax: A = laplacian(u,f,h,gx1,gxn,gy1,gyn)
%
% Inputs:
%   u: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   h: The grid spacing. Scalar.
%   gx1:    The Neumann condition on the west boundary. 2D matrix.
%   gxn:    The Neumann condition on the east boundary. 2D matrix.
%   gy1:    The Neumann condition on the south boundary. 2D matrix.
%   gyn:    The Neumann condition on the north boundary. 2D matrix.
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

    % Neumann conditions
    % West Neumann
    A(1,j,k) = (2*u(2,j,k)+u(1,j-1,k)+u(1,j+1,k)...
        +u(1,j,k-1)+u(1,j,k+1)-6*u(1,j,k)-2*h*gx1(1,j,k))./(h^2);

    % East Neumann
    A(end,j,k) = (2*u(end-1,j,k)+u(end,j-1,k)+u(end,j+1,k)...
        +u(end,j,k-1)+u(end,j,k+1)-6*u(end,j,k)+2*h*gxn(1,j,k))./(h^2);

    % South Neumann
    A(i,1,k) = (u(i-1,1,k)+u(i+1,1,k)+2*u(i,2,k)...
        +u(i,1,k-1)+u(i,1,k+1)-6*u(i,1,k)-2*h*gy1(i,1,k))./(h^2);

    % North Neumann
    A(i,end,k) = (u(i-1,end,k)+u(i+1,end,k)+2*u(i,end-1,k)...
        +u(i,end,k-1)+u(i,end,k+1)-6*u(i,end,k)+2*h*gyn(i,1,k))./(h^2);

    % South West Corner Neumann
    A(1,1,k) = (2*u(2,1,k)+2*u(1,2,k)...
        +u(1,1,k-1)+u(1,1,k+1)...
        -6*u(1,1,k)-2*h*gx1(1,1,k)-2*h*gy1(1,1,k))./(h^2);

    % South East Corner Neumann
    A(end,1,k) = (2*u(end-1,1,k)+2*u(end,2,k)...
        +u(end,1,k-1)+u(end,1,k+1)...
        -6*u(end,1,k)+2*h*gxn(1,1,k)-2*h*gy1(end,1,k))./(h^2);

    % North West Corner Neumann
    A(1,end,k) = (2*u(2,end,k)+2*u(1,end-1,k)...
        +u(1,end,k-1)+u(1,end,k+1)...
        -6*u(1,end,k)-2*h*gx1(1,end,k)+2*h*gyn(1,1,k))./(h^2);

    % North East Corner Neumann
    A(end,end,k) = (2*u(end-1,end,k)+2*u(end,end-1,k)...
        +u(end,end,k-1)+u(end,end,k+1)...
        -6*u(end,end,k)+2*h*gxn(1,end,k)+2*h*gyn(end,1,k))./(h^2);
end