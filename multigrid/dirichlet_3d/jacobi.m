function unew = jacobi(u,f,h)
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    omega = 2/3;
    N = size(u,1);
    i = 2:N-1;
    j = 2:N-1;
    k = 2:N-1;
    hsq = h^2;
    unew = u;

    % Jacobi iteration for all interior grid points
    unew(i,j,k) =(1-omega)*u(i,j,k) + omega*(1/6)*(u(i-1,j,k)+u(i+1,j,k)...
        +u(i,j-1,k)+u(i,j+1,k)+u(i,j,k-1)+u(i,j,k+1)-f(i,j,k)*hsq);
end