function unew = twogrid(u,f,nsmooth,h)
    unew = u;

    % Smooting solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h);
    end

    % Fining residual
    r = laplacian(u,h)+f;

    % Restricting residual to coarse grid
    rc = restrict(r);

    % Initial guess on error
    ec = zeros(size(rc));

    % Smoothing residual
    for i=1:nsmooth
        ec = jacobi(ec,rc,h);
    end

    % Interpolating error
    ef = interpolate(ec);

    % Updating solution
    unew = u + ef;

    % Smoothing solution
    for i=1:nsmooth
        unew = jacobi(unew,f,h);
    end
end