clear all; close all; clc;
addpath(genpath('../../../matlab_helpers'))
%% 
k = 6;
n = 2;
fig = figure('units','inch','position',[0,0,9,8]);
for l = 1:4
    pts = n^l+1;
    x = linspace(0,pi,pts);
    splt = subplot(4,1,l);
    plot(x,sin(k*x),'bo-','linewidth',2)
    ylim([-1,1])
    xlabel('$x$','interpreter','latex','fontsize',16)
    ylabel('$y$','interpreter','latex','fontsize',16)
    txt = strcat(['$',num2str(pts),'$ points']);
    title(txt,'interpreter','latex','fontsize',16)
    ax = gca;
    ax.FontSize = 10; 
end
sgtitle(strcat(['$\sin{(',num2str(k),'x)}$ Respresented on Multiple Grids']),...
    'interpreter','latex','fontsize',20)
fig2pdf(fig,'sine_on_multigrid.pdf')