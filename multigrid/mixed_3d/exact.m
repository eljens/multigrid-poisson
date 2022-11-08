function u = exact(u,f,L,U,p,h,gx1,gxn,gy1,gyn)
%EXACT
% This function computes the exact solution to \Delta u = f. This routine
% takes a precomputed LDL' factorization as input. It can be computed by
% [L,D,P,S]=ldl(A) where A is the system matrix.
%
% Syntax: u = exact(f,L,D,P,S)
%
% Inputs:
%   f: The right hand side. 3D matrix.
%   L: The lower triangular matrix from the LDL' factorization. 2D matrix.
%   D: The diagonal matrix from the LDL' factorization. 2D matrix.
%   P: Mermutation matrix from the LDL' factorization. 2D matrix.
%   S: The substitution matrix from the LDL' factorization. 2D matrix.
%
% Outputs:
%   u: The exact solution.
%
% See also: ldl
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 31st of October 2022
    N = size(f);
    f(:,:,2:end-1) = f(:,:,2:end-1)*h^2;

    % Dirichlet condition
    f(:,:,1) = u(:,:,1);
    f(:,:,end) = u(:,:,end);
    
    % West
    f(1,:,2:end-1) = f(1,:,2:end-1) + 2*h*gx1(1,:,2:end-1);

    % East
    f(end,:,2:end-1) = f(end,:,2:end-1) - 2*h*gxn(1,:,2:end-1);
    
    % South
    f(:,1,2:end-1) = f(:,1,2:end-1) + 2*h*gy1(:,1,2:end-1);
    
    % North
    f(:,end,2:end-1) = f(:,end,2:end-1) - 2*h*gyn(:,1,2:end-1);
    
    f = f(:);%reshape(f,prod(N),1);
    %u=S*P*(L'\(D\(L\((P'*S)*f))));
    %A = system_matrix(N);
    %u = A\f;
    u = U\(L\f(p));
    u = reshape(u,N);
end