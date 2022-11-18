function [X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N] = ...
    get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun)
    % number of cells in each dimension
    N = 2^l*(n-1)+1;
    
    % Axes
    x = linspace(x0,x0+len,N(1));
    y = linspace(y0,y0+len*(N(2)-1)/(N(1)-1),N(2));
    z = linspace(z0,z0+len*(N(3)-1)/(N(1)-1),N(3));

    % Creating mesh grid
    [X,Y,Z] = ndgrid(x',y',z');
    
    % The uniform grid spacing is given by
    h = x(2)-x(1);

    % Computing the true solution in order to compare results
    utrue = ufun(X,Y,Z);
    
    % Initializing the Dirichlet boaundary conditions
    u = zeros(N);
    u(1:end,1:end,1) = utrue(1:end,1:end,1);
    u(1:end,1:end,end) = utrue(1:end,1:end,end);
    
    % Intializing the right hand side
    f = ffun(X,Y,Z);
    
    % Neumann conditions
    % Western boundary
    gx1 = dxfun(X(1,:,:),Y(1,:,:),Z(1,:,:));
    
    % Eastern boundary
    gxn = dxfun(X(end,:,:),Y(end,:,:),Z(end,:,:));
    
    % Southern boundary
    gy1 = dyfun(X(:,1,:),Y(:,1,:),Z(:,1,:));
    
    % Northern boundary
    gyn = dyfun(X(:,end,:),Y(:,end,:),Z(:,end,:));
    
end