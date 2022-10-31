function A = system_matrix(n,h)
e = ones(n,n,n)/h^2;

b = e;

a = -6*e;

a(:,:,1) = 1;
a(:,:,n) = 1;
a(:,1,:) = 1;
a(:,n,:) = 1;
a(1,:,:) = 1;
a(n,:,:) = 1;
b(:,:,1) = 0;
b(:,:,n) = 0;
b(:,1,:) = 0;
b(:,n,:) = 0;
b(1,:,:) = 0;
b(n,:,:) = 0;

b = b(:);
a = a(:);

A = spdiags([b([(n*n+1):end 1:n*n]) b([(n+1):end 1:n]) b([2:end 1]) a...
    b([end 1:(end-1)]) b([(end-n+1):end 1:end-(n)]) b([(end-n*n+1):end 1:end-(n*n)])],[-n*n,-n,-1,0,1,n,n*n],n*n*n,n*n*n);

end