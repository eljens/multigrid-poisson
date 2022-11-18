function unew = Vcycle(unew,f,nsmooth,h,n,L,D,P,S)
%VCYCLE
% This function makes one Vcycle on a cubic domain with $1+2^l$ lattices in
% each dimension. Uniform grid spacing is assumed. On the coarsest domain,
% there will be only one grid point which makes it very easy to compute the
% correction.
%
% Syntax: unew = Vcycle(unew,f,nsmooth,h)
%
% Inputs:
%   unew: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   nsmooth: The number of Jacobi iterations per grid. Scalar.
%   h: The grid spacing. Scalar.
%   n: The dimensions of the coarsest domain. 1 times 3 vector.
%   L: The lower triangular matrix from the LDL' factorization. 2D matrix.
%   D: The diagonal matrix from the LDL' factorization. 2D matrix.
%   P: Mermutation matrix from the LDL' factorization. 2D matrix.
%   S: The substitution matrix from the LDL' factorization. 2D matrix.
%
% Outputs:
%   A: An approximation of the Laplacian of the solution iterate u. The
%   residuals of Au=f may be computed ad r = f-Au.
%
% See also: JACOBI, RESIDUAL, RESTRICT, INTERPOLATE, EXACT
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    assert(h > 0,"Grid spacing must be positive");
    assert(nsmooth > 0,"Number of Jacobi iterations must be at least one")
    if (size(unew,1) == n(1)) && (size(unew,2) == n(2)) && (size(unew,3) == n(3))
        unew = exact(f,L,D,P,S);
        return;
    end

    % Smooting solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h);
    end

    % Fining residual
    r = residual(unew,f,h);

    % Restricting residual to coarse grid
    rc = restrict(r);

    % Initial guess on error
    ec = zeros(size(rc));

    % Recursion
    ec = Vcycle(ec,rc,nsmooth,2*h,n,L,D,P,S);

    % Interpolating error
    ef = interpolate(ec);

    % Updating solution
    unew = unew + ef;

    % Smoothing solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h);
    end
end