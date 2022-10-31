%DRIVER
% This piece of code compares the Jacobi method with a multi-grid method
% based on Jacobi iterations on increasingly coarser grids.
%
% This implementation is simplified and will only work for cubic domains of
% size 2^l+1 where l is a positive integer.
%
% Dirichlet boundary conditions are imposed on all boundary surfaces.
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022

clear all; close all; clc;
%%
% Frequencies for test function
kx = 3;
ky = 2;
kz = 3.56;

% Test function
ufun = @(x,y,z)(sin(kx*x).*sin(ky*y).*sin(kz*z));

% Laplacian of the test function
ffun = @(x,y,z)(-(kx^2+ky^2+kz^2)*sin(kx*x).*sin(ky*y).*sin(kz*z));

% Domian offsets from origo
x0 = -0.4;
y0 = 0.21;
z0 = 1.23;

% Length of domain sides
len = pi;

% Number of layers in the multi-grid
l = 6;

% number of cells in each dimension
N = 2^l+1;

% Number of smoothings per grids
nsmooth = 15;

% Maximal number of iterations
max_iter = 4000;
max_time = 150; % seconds

% Tolerance
tol = 100*eps;

% Axes
x = linspace(x0,x0+len,N);
y = linspace(y0,y0+len,N);
z = linspace(z0,z0+len,N);

% Creating mesh grid
[X,Y,Z] = meshgrid(x',y',z');

% The uniform grid spacing is given by
h = x(2)-x(1);

% Initializing the Dirichlet boaundary conditions
u = ufun(X,Y,Z);
u(2:end-1,2:end-1,2:end-1) = zeros(N-2,N-2,N-2);

% Intializing the right hand side
f = ffun(X,Y,Z);

% Computing the true solution in order to compare results
utrue = ufun(X,Y,Z);

%% Jacobi solver

% Copying solution array
u_jac = zeros(size(u))+u;

% Vector for storing residuals
err_jac = zeros(max_iter,1);
t_jac = zeros(max_iter,1);

start = tic;
for i=1:max_iter
    u_jac = jacobi(u_jac,f,h);
    r = residual(u_jac,f,h);
    err_jac(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);
    t_jac(i) = toc(start);
    if (err_jac(i) < tol) || (t_jac(i) > max_time)
        err_jac = err_jac(1:i);
        t_jac = t_jac(1:i);
        break;
    end
end

%% Multigrid solver
err_vcycle = zeros(max_iter,1);
u_vcycle = zeros(size(u))+u;
t_vcycle = zeros(max_iter,1);

start = tic;
for i=1:max_iter
    u_vcycle = Vcycle(u_vcycle,f,nsmooth,h);
    r = -residual(u_vcycle,f,h);
    err_vcycle(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);%max(max(abs(u_vcycle-utrue)));
    t_vcycle(i) = toc(start);
    if (err_vcycle(i) < tol) || (t_vcycle(i) > max_time)
        err_vcycle = err_vcycle(1:i);
        t_vcycle = t_vcycle(1:i);
        break;
    end
end

%% Plot of some slice of the solution

plotidx = floor(N/2);

fig = figure('units','inch','position',[0,0,15,4]);
subplot(1,3,1)
surf(X(:,:,plotidx),Y(:,:,plotidx),utrue(:,:,plotidx))
title("True Solution")
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
view([0 0 90])
colorbar
subplot(1,3,2)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_jac(:,:,plotidx))
title("Jacobi Approximation")
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
view([0 0 90])
colorbar
subplot(1,3,3)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_vcycle(:,:,plotidx))
title("Vcycle Approximation")
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