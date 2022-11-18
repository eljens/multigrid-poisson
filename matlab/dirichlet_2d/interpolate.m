function uf = interpolate(uc)
    n = size(uc,1);
    uf = zeros(2*n-1,2*n-1);
    uf(1:2:end,1:2:end) = uc;
    uf(2:2:end,:) = (uf(3:2:end,:) + uf(1:2:end-2,:))/2;
    uf(:,2:2:end) = (uf(:,3:2:end) + uf(:,1:2:end-2))/2;
end