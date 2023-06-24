clear all; close all; clc;

%% 
gpu = 'V100';%'MI250X';

f1 = figure('units','inch','position',[0,0,6,4]);
for l=1:2:7
    data = parseFile(strcat([gpu,'/weak_scaling_large_',num2str(l),'.txt']));
    avgtime = data.seconds./data.maxiter;
    plot(data.devices,avgtime(1)./avgtime,'*--','DisplayName',...
    strcat(['$l=',num2str(l),'$']),'linewidth',2)
    if l==1
        hold on
    end
end
hold off
grid()
legend('interpreter','latex','fontsize',14,'location','southwest')
xlabel('Number of Devices','interpreter','latex','fontsize',20)
ylabel('Efficiency','interpreter','latex','fontsize',20)
if strcmp(gpu,'V100')
    ylim([0.8 1.05])
    xticks([1,2,4])
else 
    ylim([0.6 1.05])
    xticks([1,2,4,8])
end
set(gca,'FontSize',16)
saveas(gcf,strcat(['../figures/weak_scaling_',gpu,'.png']))
fig2pdf(gca,strcat(['weak_scaling_',gpu,'.pdf']))

%% Reading content in file
function data = parseFile(filename)
    formatSpec = '%f %f %f %f %f %f %f %f %f %f';
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