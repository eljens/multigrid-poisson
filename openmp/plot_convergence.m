clear all; close all; clc;
addpath(genpath('../matlab/mixed_3d'))
%% Reading content in file
formatSpec = '%f %f %f %d %f %d %d %d';
% The eight variables are
% threads devices version N ite tol wall_time max_err

fid = fopen('results/vcycle.txt','r');
A = textscan(fid,formatSpec,'CommentStyle','#');
fclose(fid);

%% Reading variable names
fid = fopen('results/vcycle.txt','r');
varnames = strsplit(fgetl(fid), ' ');
fclose(fid);
% Deleting comment style
varnames(1) = [];

%% Accumulating data
data = struct;
for i=1:size(A,2)
    data.(varnames{i}) = A{i};
end

%% 

[alpha,beta] = ols_log_fit(data.abs_err(2:end),data.spacing(2:end));

figure(1)
loglog(data.spacing,data.abs_err,'b*-','DisplayName',...
    strcat(['Abs Err $O(h^{',num2str(beta),'})$']),'linewidth',2)
hold on
loglog(data.spacing,8*data.spacing.^2,'k--',...
    'DisplayName','$\mathcal{O}(h^2)$','linewidth',2)
hold off
grid()
legend('interpreter','latex','fontsize',14,'location','nw')
xlabel('$h$','interpreter','latex','fontsize',18)
ylabel('Absolute Error','interpreter','latex','fontsize',18)
saveas(gcf,'./figures/mg_convergence.png')
