function [X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N] = ...
    get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun)
%GET_DOMAIN
% This function returns both the grid and the boundary conditions together
% with initialized matrices for the solution iterate, the right hand side,
% the true solution, the grid spacing and the grid dimensions.
%
% Syntax: [X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N]
%         = get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun)
%
% Inputs:
%   n:      The grid dimensions of the coarsest grid. Scalar
%   l:      The number of grids in the Vcycle. Scalar
%   len:    The Length of the x axis. Scalar.
%   x0:     The x component of the origin. Scalar.
%   y0:     The y component of the origin. Scalar.
%   z0:     The z component of the origin. Scalar.
%   ufun:   Function to compute the true solution.
%   ffun:   Function to compute the right hand side.
%   dxfun:  Function to compute the x derivative of u(x,y,z).
%   dyfun:  Function to compute the y derivative of u(x,y,z).
%
% Outputs:
%   X:      The x component of the meshgrid. 3d matrix.
%   Y:      The y component of the meshgrid. 3d matrix.
%   Z:      The z component of the meshgrid. 3d matrix.
%   gx1:    The Neumann condition on the west boundary. 2D matrix.
%   gxn:    The Neumann condition on the east boundary. 2D matrix.
%   gy1:    The Neumann condition on the south boundary. 2D matrix.
%   gyn:    The Neumann condition on the north boundary. 2D matrix.
%   u:      The initial solution iterate with Dirichlet conditions imposed. 
%           3D matrix.
%   f:      The right-hand side. 3D matrix.
%   utrue:  The true solution to the test problem. 3D matrix.
%   h:      The uniform grid spacing. Scalar.
%   N:      The grid dimensions for the finest grid. 1 times 3 vector.
%
% See also: problem_definition
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 18th of November 2022
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