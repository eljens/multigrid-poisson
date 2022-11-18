function [u,t,abse,relres] = multigrid_solver(u,f,utrue,nsmooth,h,n,gx1,gxn,gy1,gyn,max_iter)
    
    % Precomputing LU  factorization for coarsest grid
    A = system_matrix(n);
    [L,U,p] = lu(A,'vector');

    abse = zeros(max_iter,1);
    relres = zeros(max_iter,1);
    t = zeros(max_iter,1);
    
    start = tic;
    for i=1:max_iter
        u = Vcycle(u,f,nsmooth,h,n,L,U,p,gx1,gxn,gy1,gyn);
        r = residual(u,f,h,gx1,gxn,gy1,gyn);
        t(i) = toc(start);
        relres(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);
        abse(i) = max(max(max(abs(u-utrue))));
    end
end