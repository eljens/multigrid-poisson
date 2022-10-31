clear all; close all; clc;
%%
kx = 3;
ky = 2;
kz = 3.56;
ufun = @(x,y,z)(sin(kx*x).*sin(ky*y).*sin(kz*z));
ffun = @(x,y,z)(-(kx^2+ky^2+kz^2)*sin(kx*x).*sin(ky*y).*sin(kz*z));

% Domian size
x0 = -0.4;
y0 = 0.21;
z0 = 1.23;
len = pi;

% 
l = 5;

% number of cells
N = 2^l+1;

% Number of smoothings
nsmooth = 15;

% Maximal number of iterations
max_iter = 1000;

% Axes
x = linspace(x0,x0+len,N);
y = linspace(y0,y0+len,N);
z = linspace(z0,z0+len,N);

% grid
[X,Y,Z] = meshgrid(x',y',z');

% grid spacing
h = x(2)-x(1);

% Domain
u = ufun(X,Y,Z);
u(2:end-1,2:end-1,2:end-1) = zeros(N-2,N-2,N-2);
f = ffun(X,Y,Z);

% True solution
utrue = ufun(X,Y,Z);

%% Jacobi solver

u_jac = zeros(size(u))+u;
err_jac = zeros(max_iter,1);
for i=1:max_iter
    u_jac = jacobi(u_jac,f,h);
    r = -laplacian(u_jac,f,h)+f;
    err_jac(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);%max(max(abs(u_jac-utrue)));
end

%% Multigrid solver
err_vcycle = zeros(max_iter,1);
u_vcycle = zeros(size(u))+u;

for i=1:max_iter
    u_vcycle = Vcycle(u_vcycle,f,nsmooth,h);
    r = -laplacian(u_vcycle,f,h)+f;
    err_vcycle(i) = norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2);%max(max(abs(u_vcycle-utrue)));
    if (err_vcycle(i) < 100*eps)
        err_vcycle = err_vcycle(1:i);
        break;
    end
end

%%

plotidx = floor(N/2);

figure(2)
subplot(1,3,1)
surf(X(:,:,plotidx),Y(:,:,plotidx),utrue(:,:,plotidx))
subplot(1,3,2)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_jac(:,:,plotidx))
subplot(1,3,3)
surf(X(:,:,plotidx),Y(:,:,plotidx),u_vcycle(:,:,plotidx))

%%

figure(3)
semilogy(1:max_iter,err_jac,'b.','displayname','Jacobi')
hold on
semilogy(1:size(err_vcycle,1),err_vcycle,'r.','displayname','Vcycle')
hold off
legend('interpreter','latex','fontsize',16,'location','east')
xlabel('Iterations','interpreter','latex','fontsize',18)
ylabel('$\frac{||f-Au||_2}{||f||_2}$','interpreter','latex','fontsize',18)
grid()
title('Convergence Study','interpreter','latex','fontsize',20)