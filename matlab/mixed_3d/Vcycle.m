function unew = Vcycle(unew,f,nsmooth,h,n,L,U,p,gx1,gxn,gy1,gyn)
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
        unew = exact(unew,f,L,U,p,h,gx1,gxn,gy1,gyn);
        return;
    end

    % Smooting solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h,gx1,gxn,gy1,gyn);
    end

    % Fining residual
    r = residual(unew,f,h,gx1,gxn,gy1,gyn);

    % Restricting residual to coarse grid
    rc = restrict(r);

    % Initial guess on error
    ec = zeros(size(rc));
    %ec(:,:,1) = unew(1:2:end,1:2:end,1);
    %ec(:,:,end) = unew(1:2:end,1:2:end,end);

    % Updating boundary conditions
%     gx1c = gx1 - (unew(2,:,:)-unew(1,:,:))./h;
%     gxnc = gxn - (unew(end,:,:)-unew(end-1,:,:))./h;
%     gy1c = gy1 - (unew(:,2,:)-unew(:,1,:))./h;
%     gync = gyn - (unew(:,end,:)-unew(:,end-1,:))./h;
    gx1c = gx1 - (-0.5*unew(3,:,:)+2*unew(2,:,:)-(3/2)*unew(1,:,:))./h;
    gxnc = gxn - ((3/2)*unew(end,:,:)-2*unew(end-1,:,:)+0.5*unew(end-2,:,:))./h;
    gy1c = gy1 - (-0.5*unew(:,3,:)+2*unew(:,2,:)-(3/2)*unew(:,1,:))./h;
    gync = gyn - ((3/2)*unew(:,end,:)-2*unew(:,end-1,:)+0.5*unew(:,end-2,:))./h;
    gx1c = gx1c(:,1:2:end,1:2:end);
    gxnc = gxnc(:,1:2:end,1:2:end);
    gy1c = gy1c(1:2:end,:,1:2:end);
    gync = gync(1:2:end,:,1:2:end);

    % Recursion
    ec = Vcycle(ec,rc,nsmooth,2*h,n,L,U,p,...
        gx1c,gxnc,gy1c,gync);

    % Interpolating error
    ef = interpolate(ec);

    % Updating solution
    %unew(:,:,2:end-1) = unew(:,:,2:end-1) + ef(:,:,2:end-1);
    unew = unew + ef;

    % Smoothing solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h,gx1,gxn,gy1,gyn);
    end
end