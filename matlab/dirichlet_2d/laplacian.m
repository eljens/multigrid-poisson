function A = laplacian(u,f,h)
    N = size(u,1);
    i = 2:N-1;
    j = 2:N-1;
    A = f;
    A(i,j) = (u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-4*u(i,j))./(h^2);
end