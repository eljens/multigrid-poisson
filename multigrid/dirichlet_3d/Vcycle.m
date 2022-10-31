function unew = Vcycle(unew,f,nsmooth,h)
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
%
% Outputs:
%   A: An approximation of the Laplacian of the solution iterate u. The
%   residuals of Au=f may be computed ad r = f-Au.
%
% See also: JACOBI, RESIDUAL, RESTRICT, INTERPOLATE
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    assert(h > 0,"Grid spacing must be positive");
    assert(nsmooth > 0,"Number of Jacobi iterations must be at least one")
    if(size(unew,1) == 3)
        unew = f(2,2,2);
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
    ec = Vcycle(ec,rc,nsmooth,h);

    % Interpolating error
    ef = interpolate(ec);

    % Updating solution
    unew = unew + ef;

    % Smoothing solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h);
    end
end