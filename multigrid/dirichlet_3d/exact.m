function u = exact(f,h,L,D,P,S)
    n = size(f,1);
    A = system_matrix(n,h);
    f = reshape(f,[],1);
    u=S*P*(L'\(D\(L\((P'*S)*f))));
    u = reshape(u,[n,n,n]);
end