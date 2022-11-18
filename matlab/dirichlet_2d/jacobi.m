function unew = jacobi(u,f,h)
    omega = 2/3;
    N = size(u,1);
    i = 2:N-1;
    j = 2:N-1;
    hsq = h^2;
    unew = u;
    unew(i,j) =(1-omega)*u(i,j) + omega*0.25*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-f(i,j)*hsq);
end