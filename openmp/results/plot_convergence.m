clear all; close all; clc;

%% 
problems = {'neumann','dirichlet'};
for i = 1:2
    prob = problems{i};
prob1 = parseFile(strcat(['MI250X/',prob,'_1.txt']));
prob2 = parseFile(strcat(['MI250X/',prob,'_2.txt']));
prob3 = parseFile(strcat(['MI250X/',prob,'_3.txt']));

[alpha1,beta1] = ols_log_fit(prob1.abs_err(3:end),prob1.spacing(3:end));
[alpha2,beta2] = ols_log_fit(prob2.abs_err(3:end),prob2.spacing(3:end));
[alpha3,beta3] = ols_log_fit(prob3.abs_err(3:end),prob3.spacing(3:end));

f1 = figure('units','inch','position',[0,0,6,4]);
loglog(prob1.spacing(2:end),prob1.abs_err(2:end),'*--','DisplayName',...
    strcat(['Problem 1 $ROC=',num2str(beta1),'$']),'linewidth',3,'markersize',8)
hold on
%loglog(prob2.spacing,prob2.abs_err,'^--','DisplayName',...
%    strcat(['Problem 2 $ROC=',num2str(beta2),'$']),'linewidth',2)
loglog(prob3.spacing(2:end),prob3.abs_err(2:end),'x--','DisplayName',...
    strcat(['Problem 3 $ROC=',num2str(beta3),'$']),'linewidth',3,'markersize',10)
loglog(prob1.spacing(2:end),10*prob1.spacing(2:end).^2,'k--',...
    'DisplayName','$\mathcal{O}(h^2)$','linewidth',3)
hold off
grid()
legend('interpreter','latex','fontsize',16,'location','nw')
xlabel('$h$','interpreter','latex','fontsize',22)
ylabel('Absolute Error','interpreter','latex','fontsize',22)
ylim([1e-10 5])
set(gca,'FontSize',16)
saveas(gcf,'../figures/mg_convergence.png')
fig2pdf(gca,strcat([prob,'_convergence.pdf']))
end

%% Reading content in file
function data = parseFile(filename)
    formatSpec = '%f %f %f %d %f %d %d %d %d %d';
    % The eight variables are
    % threads devices version N ite tol wall_time max_err
    
    fid = fopen(filename,'r');
    A = textscan(fid,formatSpec,'CommentStyle','#');
    fclose(fid);
    
    %% Reading variable names
    fid = fopen(filename,'r');
    varnames = strsplit(fgetl(fid), ' ');
    fclose(fid);
    % Deleting comment style
    varnames(1) = [];
    
    %% Accumulating data
    data = struct;
    for i=1:size(A,2)
        data.(varnames{i}) = A{i};
    end
end