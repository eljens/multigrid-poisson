function unew = Vcycle(unew,f,nsmooth,h,n,L,U,p,gx1,gxn,gy1,gyn)
%VCYCLE
% This function makes one Vcycle on a hyperrectangular domain with 
% 1+2^l*n lattices in each dimension. Uniform grid spacing is assumed. The
% LU decomposition of the Laplacian operator on the coarsest grid can be
% computed with [L,U,p] = lu(A,'vector').
%
% Syntax: unew = Vcycle(unew,f,nsmooth,h,n,L,U,p,gx1,gxn,gy1,gyn)
%
% Inputs:
%   unew: The current iterate of the solution. 3D matrix.
%   f: The right hand side. 3D matrix.
%   nsmooth: The number of Jacobi iterations per grid. Scalar.
%   h: The grid spacing. Scalar.
%   n: The dimensions of the coarsest domain. 1 times 3 vector.
%   L: The lower triangular matrix from the LU decomposition. 2D matrix.
%   D: The upper triangular matrix from the U decomposition. 2D matrix.
%   p: Mermutation matrix from the LU decomposition. 1D vector.
%   gx1: The Neumann condition on the west boundary. 2D matrix.
%   gxn: The Neumann condition on the east boundary. 2D matrix.
%   gy1: The Neumann condition on the south boundary. 2D matrix.
%   gyn: The Neumann condition on the north boundary. 2D matrix.
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
    % If we have reached the coarsest level, we solve Ae=r using a direct
    % method, in this case based on the LU decomposition of the discrete
    % Laplacian operator
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

    % Updating boundary conditions
    % A three point one-sided approximation of the first order derivative
    % is used to find the defect in the Neumann condition which will be the
    % Neumann condition on the coarser grid
    gx1c = gx1 - (-0.5*unew(3,:,:)+2*unew(2,:,:)-(3/2)*unew(1,:,:))./h;
    gxnc = gxn - ((3/2)*unew(end,:,:)-2*unew(end-1,:,:)+0.5*unew(end-2,:,:))./h;
    gy1c = gy1 - (-0.5*unew(:,3,:)+2*unew(:,2,:)-(3/2)*unew(:,1,:))./h;
    gync = gyn - ((3/2)*unew(:,end,:)-2*unew(:,end-1,:)+0.5*unew(:,end-2,:))./h;
    gx1c = gx1c(:,1:2:end,1:2:end);
    gxnc = gxnc(:,1:2:end,1:2:end);
    gy1c = gy1c(1:2:end,:,1:2:end);
    gync = gync(1:2:end,:,1:2:end);

    % Recursion
    % On the next grid the grid spacing will be 2h
    % On the next grid we solve Ae = r where e is the additive error and r
    % is the residual/defect. On the finest grid we solver Au=f
    ec = Vcycle(ec,rc,nsmooth,2*h,n,L,U,p,...
        gx1c,gxnc,gy1c,gync);

    % Interpolating the error from the coarser to the finer grid with
    % trilinear interpolation
    ef = interpolate(ec);

    % Updating solution by adding the error to the solution iterate
    unew = unew + ef;

    % Smoothing the high frequent errors introduced on the fine grid after
    % adding the interpolated error from the coarse grid
    for i=1:nsmooth
        unew = jacobi(unew,f,h,gx1,gxn,gy1,gyn);
    end
end