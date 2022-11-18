clear all; close all; clc;
%% TESTIING CONVERGENCE OF ENTIRE VCYCLE FOR DECREASING GRID SPACING

[ufun,ffun,dxfun,dyfun] = problem_definition();

% Domain offsets from origo
x0 = -0.71;
y0 = 0.21;
z0 = 1.23;

% Length of domain sides
len = pi/20;

% Number of smoothings per grids
nsmooth = 20;

% Maximal number of iterations
max_iter = 10;
max_time = 1800; % seconds

% Tolerance
tol = 100*eps;

% Number of cells in coarsest grid - must be odd
n = [7,3,3];
assert(sum(mod(n,2) == 1) == 3,"Number of elements on coarsest grid must be odd")

%% Looping over l

err_vec = [];
h_vec = [];

for l=3:6
    [X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N] = get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun);

    % Calling multigrid solver
    [u,t,abse,relres] = multigrid_solver(u,f,utrue,nsmooth,h,n,gx1,gxn,gy1,gyn,max_iter);

    err_vec = [err_vec;abse(end)];
    h_vec = [h_vec;h];
    disp(strcat(['Solved on level ',num2str(l)]))
end
%%
[alpha,beta] = ols_log_fit(err_vec(2:end),h_vec(2:end));

figure(1)
loglog(h_vec,err_vec,'b*-','DisplayName',...
    strcat(['Abs Err $O(h^{',num2str(beta),'})$']),'linewidth',2)
hold on
loglog(h_vec,12*h_vec.^2,'k--',...
    'DisplayName','$\mathcal{O}(h^2)$','linewidth',2)
hold off
grid()
legend('interpreter','latex','fontsize',14,'location','nw')
xlabel('$h$','interpreter','latex','fontsize',18)
ylabel('Absolute Error','interpreter','latex','fontsize',18)
saveas(gcf,'./figures/mg_convergence.png')