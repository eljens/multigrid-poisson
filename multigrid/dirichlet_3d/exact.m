function u = exact(f,L,D,P,S)
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
    f = reshape(f,[],1);
    u=S*P*(L'\(D\(L\((P'*S)*f))));
    u = reshape(u,N);
end