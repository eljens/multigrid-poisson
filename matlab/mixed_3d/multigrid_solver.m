function [u,t,abse,relres] = ...
    multigrid_solver(u,f,utrue,nsmooth,h,n,gx1,gxn,gy1,gyn,maxiter)
%MULTIGRID_SOLVER
% This function runs maxiter iterations of the Vcycle and returns the
% result of the simulation together with absolute error and relative
% tolerance measures.
%
% Syntax: [u,t,abse,relres] = 
%         multigrid_solver(u,f,utrue,nsmooth,h,n,gx1,gxn,gy1,gyn,maxiter)
%
% Inputs:
%   u:      The current iterate of the solution. 3D matrix.
%   f:      The right hand side. 3D matrix.
%   utrue   The true solution. Used to measure the absolute error. 3D
%           matrix.
%   nsmooth:The number of Jacobi iterations per grid. Scalar.
%   h:      The grid spacing. Scalar.
%   n:      The dimensions of the coarsest grid. 1 times 3 vector.
%   gx1:    The Neumann condition on the west boundary. 2D matrix.
%   gxn:    The Neumann condition on the east boundary. 2D matrix.
%   gy1:    The Neumann condition on the south boundary. 2D matrix.
%   gyn:    The Neumann condition on the north boundary. 2D matrix.
%   maxiter:The maximum number of iterations. Scalar.
%
% Outputs:
%   u:      The resulting approximate solution. 3D matrix.
%   t:      The time taken for each iteration. maxiter times 1 vector.
%   abse:   The absolute error of the approximate solution. 
%           maxiter times 1 vector.
%   relres: The relative residuals. maxiter times 1 vector.
%
% See also: Vcycle, LU, residual, system_matrix
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 18th of November 2022
    % Precomputing LU  factorization for coarsest grid
    A = system_matrix(n);
    [L,U,p] = lu(A,'vector');

    abse = zeros(maxiter,1);
    relres = zeros(maxiter,1);
    t = zeros(maxiter,1);
    
    start = tic;
    for i=1:maxiter
        % Running a Vcycle iteration
        u = Vcycle(u,f,nsmooth,h,n,L,U,p,gx1,gxn,gy1,gyn);

        % Evaluating the solution and storing the relative residuals and
        % the maximal error measured in the infinity norm.
        r = residual(u,f,h,gx1,gxn,gy1,gyn);
        t(i) = toc(start);
        relres(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);
        abse(i) = max(max(max(abs(u-utrue))));
    end
end