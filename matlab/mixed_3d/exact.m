function u = exact(u,f,L,U,p,h,gx1,gxn,gy1,gyn)
%EXACT
% This function computes the exact solution to \Delta u = f. This routine
% takes a precomputed LU decomposition as input. The LU decomposition of 
% the Laplacian operator on the coarsest grid can be computed with 
% [L,U,p] = lu(A,'vector').
%
% Syntax: u = exact(u,f,L,U,p,h,gx1,gxn,gy1,gyn)
%
% Inputs:
%   u: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   L: The lower triangular matrix from the LU decomposition. 2D matrix.
%   D: The upper triangular matrix from the U decomposition. 2D matrix.
%   p: Mermutation matrix from the LU decomposition. 1D vector.
%   h: The uniform grid spacing. Scalar.
%   gx1: The Neumann condition on the west boundary. 2D matrix.
%   gxn: The Neumann condition on the east boundary. 2D matrix.
%   gy1: The Neumann condition on the south boundary. 2D matrix.
%   gyn: The Neumann condition on the north boundary. 2D matrix.
%
% Outputs:
%   u: The exact solution.
%
% See also: LU
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 31st of October 2022
    N = size(f);
    f(:,:,2:end-1) = f(:,:,2:end-1)*h^2;

    % Imposing the Dirichlet conditions
    % Bottom
    f(:,:,1) = u(:,:,1);

    % Top
    f(:,:,end) = u(:,:,end);
    
    % Imposing the Neumann conditions
    % West
    f(1,:,2:end-1) = f(1,:,2:end-1) + 2*h*gx1(1,:,2:end-1);

    % East
    f(end,:,2:end-1) = f(end,:,2:end-1) - 2*h*gxn(1,:,2:end-1);
    
    % South
    f(:,1,2:end-1) = f(:,1,2:end-1) + 2*h*gy1(:,1,2:end-1);
    
    % North
    f(:,end,2:end-1) = f(:,end,2:end-1) - 2*h*gyn(:,1,2:end-1);
    
    % Flattening the right hand side
    f = f(:);
    
    % Solving for the exact solution using the precomputed LU factorization
    u = U\(L\f(p));

    % Reshaping the solution from1D to 3D
    u = reshape(u,N);
end