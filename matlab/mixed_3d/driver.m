%DRIVER
% This piece of code compares the Jacobi method with a multi-grid method
% based on Jacobi iterations on increasingly coarser grids.
%
% This implementation assumes that the finest grid will always have an even
% number of cells, or equivalently an odd number of grid points in each
% dimension. This is only for easing the implementation of the
% interpolation and restriction operations.
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022

clear all; close all; clc;
%%
% Test problem
[ufun,ffun,dxfun,dyfun] = problem_definition();

% Domain offsets from origo
x0 = -0.71;
y0 = 0.21;
z0 = 1.23;

% Length of domain sides
len = pi+0.5;

% Number of layers in the multi-grid
l = 5;

% Number of cells in coarsest grid - must be odd
n = [5, 5, 5];
assert(sum(mod(n,2) == 1) == 3,"Number of elements on coarsest grid must be odd")

% Number of smoothings per grids
nsmooth = 20;

% Maximal number of iterations
max_iter = 2000;
max_time = 20; % seconds

% Tolerance
tol = 100*eps;

% Initializing the domain and the test problem
[X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N] = ...
    get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun);

%%
if prod(N) < 25000
    A = system_matrix(N);
    [L,U,p] = lu(A,'vector');
    usolve = exact(u,f,L,U,p,h,gx1,gxn,gy1,gyn);
    
    idx = round(N(3)/2);
    figure(1)
    subplot(1,3,1)
    surf(X(:,:,idx),Y(:,:,idx),utrue(:,:,idx))
    xlabel('x')
    ylabel('y')
    subplot(1,3,2)
    surf(X(:,:,idx),Y(:,:,idx),usolve(:,:,idx))
    xlabel('x')
    ylabel('y')
    subplot(1,3,3)
    surf(X(:,:,idx),Y(:,:,idx),utrue(:,:,idx)-usolve(:,:,idx))
    xlabel('x')
    ylabel('y')
end

%% Jacobi solver

% Copying solution array
u_jac = zeros(size(u))+u;

% Vector for storing residuals
err_jac = zeros(max_iter,1);
t_jac = zeros(max_iter,1);

start = tic;
for i=1:max_iter
    u_jac = jacobi(u_jac,f,h,gx1,gxn,gy1,gyn);
    r = residual(u_jac,f,h,gx1,gxn,gy1,gyn);
    err_jac(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);
    t_jac(i) = toc(start);
    if (err_jac(i) < tol) || (t_jac(i) > max_time)
        err_jac = err_jac(1:i);
        t_jac = t_jac(1:i);
        break;
    end
end

%% Multigrid solver
u_vcycle = zeros(size(u))+u;
maxiter = 10;

% Calling the Vcycle solver which uses a LU factorization on the coarsest
% level
[u_vcycle,t_vcycle,err_vcycle,relres] = ...
    multigrid_solver(u_vcycle,f,utrue,nsmooth,h,n,gx1,gxn,gy1,gyn,maxiter);

%% Plot of some slice of the solution

plotidx = floor(N(3)/2);

fig = figure('units','inch','position',[0,0,15,15]);
subplot(2,2,1)
surf(X(:,:,plotidx),Y(:,:,plotidx),utrue(:,:,plotidx))
title("True Solution")
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
view([0 0 90])
colorbar
% Exact solution can only be computed for small problems in reasonable time
if prod(N) < 25000
    subplot(2,2,2)
    surf(X(:,:,plotidx),Y(:,:,plotidx),usolve(:,:,plotidx))
    title("LU Solution")
    xlabel('$x$','interpreter','latex','fontsize',16)
    ylabel('$y$','interpreter','latex','fontsize',16)
    view([0 0 90])
    colorbar
end
subplot(2,2,3)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_jac(:,:,plotidx))
title("Jacobi Approximation")
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
view([0 0 90])
colorbar
subplot(2,2,4)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_vcycle(:,:,plotidx))
title("Vcycle Approximation")
%surf(X(:,:,plotidx),Y(:,:,plotidx),utrue(:,:,plotidx)-u_jac(:,:,plotidx))
%title("Error")
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
view([0 0 90])
colorbar
saveas(fig,"./figures/solution.png")

%% Plot of tolerance as function of iteration

figure(3)
semilogy(1:size(err_jac,1),err_jac,'b.','displayname','Jacobi')
hold on
semilogy(1:size(err_vcycle,1),err_vcycle,'r.','displayname','Vcycle')
hold off
legend('interpreter','latex','fontsize',16,'location','northeast')
xlabel('Iterations','interpreter','latex','fontsize',18)
ylabel('$\frac{||f-Au||_2}{||f||_2}$','interpreter','latex','fontsize',18)
grid()
title('Convergence Study: Iterations','interpreter','latex','fontsize',20)
saveas(gcf,'./figures/iterations.png')

%% Plot of tolerance as function of time

max_t = min(t_jac(end),t_vcycle(end));

figure(4)
semilogy(t_jac(t_jac < max_t),err_jac(t_jac < max_t),'b.',...
    'displayname','Jacobi')
hold on
semilogy(t_vcycle(t_vcycle < max_t),err_vcycle(t_vcycle < max_t),...
    'r.','displayname','Vcycle')
hold off
legend('interpreter','latex','fontsize',16,'location','northeast')
xlabel('Time in Seconds','interpreter','latex','fontsize',18)
ylabel('$\frac{||f-Au||_2}{||f||_2}$','interpreter','latex','fontsize',18)
grid()
title('Convergence Study: Time','interpreter','latex','fontsize',20)
saveas(gcf,'./figures/time.png')