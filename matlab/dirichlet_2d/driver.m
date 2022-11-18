clear all; close all; clc;
%%
kx = 3;
ky = 2;
ufun = @(x,y)(sin(kx*x).*sin(ky*y));
ffun = @(x,y)(-(kx^2+ky^2)*sin(kx*x).*sin(ky*y));

% Domian size
x0 = -0.4;
y0 = 0.21;
len = pi;

% 
l = 5;

% number of cells
N = 2^l+1;

% Number of smoothings
nsmooth = 10;

% Maximal number of iterations
max_iter = 500;

% Axes
x = linspace(x0,x0+len,N);
y = linspace(y0,y0+len,N);

% grid
[X,Y] = meshgrid(x',y');

% grid spacing
h = x(2)-x(1);

% Domain
u = ufun(X,Y);
u(2:end-1,2:end-1) = zeros(N-2,N-2);
f = ffun(X,Y);

% True solution
utrue = ufun(X,Y);

%% Jacobi Solver

u_jac = zeros(size(u))+u;
err_jac = zeros(max_iter,1);
for i=1:max_iter
    u_jac = jacobi(u_jac,f,h);
    r = -laplacian(u_jac,f,h)+f;
    err_jac(i) = norm(r,2)/norm(f,2);%max(max(abs(u_jac-utrue)));
end

%% Multi-Grid Solver

err_vcycle = zeros(max_iter,1);
u_vcycle = zeros(size(u))+u;

for i=1:max_iter
    u_vcycle = Vcycle(u_vcycle,f,nsmooth,h);
    r = -laplacian(u_vcycle,f,h)+f;
    err_vcycle(i) = norm(r,2)/norm(f,2);%max(max(abs(u_vcycle-utrue)));
    if (0)
        figure(4)
        surf(X,Y,utrue-u_vcycle)
        pause(1)
    end
end

%%
figure(2)
subplot(1,3,1)
surf(X,Y,utrue)
subplot(1,3,2)
surf(X,Y,u_jac)
subplot(1,3,3)
surf(X,Y,u_vcycle)

%% 
figure(3)
semilogy(1:max_iter,err_jac,'displayname','Jacobi')
hold on
semilogy(1:max_iter,err_vcycle,'displayname','Vcycle')
hold off
legend()
