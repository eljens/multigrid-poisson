clear all; close all; clc;
%% TESTIING CONVERGENCE OF ENTIRE VCYCLE FOR DECREASING GRID SPACING

[ufun,ffun,dxfun,dyfun] = problem_definition();

% Domain offsets from origo
x0 = -0.71;
y0 = 0.21;
z0 = 1.23;

% Length of domain sides
len = pi/10;

% Number of smoothings per grids
nsmooth = 20;

% Maximal number of iterations
max_iter = 10;
max_time = 1800; % seconds

% Tolerance
tol = 100*eps;

% Number of cells in coarsest grid - must be odd
n = [3,3,5];
assert(sum(mod(n,2) == 1) == 3,"Number of elements on coarsest grid must be odd")

%% Looping over l

abse = [];
relres = [];
h_vec = [];

for l=1:3
    [X,Y,Z,gx1,gxn,gy1,gyn,u,f,utrue,h,N] = get_domain(n,l,len,x0,y0,z0,ufun,ffun,dxfun,dyfun);
    A = system_matrix(N);
    %A(3,3) = 3;
    [L,U,p] = lu(A,'vector');
    u = exact(u,f,L,U,p,h,gx1,gxn,gy1,gyn);
    r = residual(u,f,h,gx1,gxn,gy1,gyn);
    relres = [relres;norm(reshape(r,[],1),2)/norm(reshape(f,[],1),2)];
    abse = [abse;max(max(max(abs(u-utrue))))];
    h_vec = [h_vec;h];
    disp(strcat(['Finished iteration ',num2str(l)]))
end
%%
[alpha,beta] = ols_log_fit(abse,h_vec);

figure(1)
%loglog(h_vec,relres,'DisplayName','Relative Residual')
%hold on
loglog(h_vec,abse,'b*-','DisplayName',...
    strcat(['Abs Err $O(h^{',num2str(beta),'})$']),'linewidth',2)
hold on
loglog(h_vec,h_vec.^2,'k--','DisplayName',...
   '$\mathcal{O}(h^2)$','linewidth',2)
hold off
grid()
legend('interpreter','latex','fontsize',14,'location','nw')
xlabel('$h$','interpreter','latex','fontsize',18)
ylabel('Absolute Error','interpreter','latex','fontsize',18)
saveas(gcf,'./figures/sys_mat_convergence.png')