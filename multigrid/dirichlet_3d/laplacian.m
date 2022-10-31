function A = laplacian(u,f,h)
    N = size(u,1);
    i = 2:N-1;
    j = 2:N-1;
    k = 2:N-1;
    A = f;
    A(i,j,k) = (u(i-1,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i,j+1,k)...
        +u(i,j,k-1)+u(i,j,k+1)-6*u(i,j,k))./(h^2);
end