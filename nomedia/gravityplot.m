clear all
close all

 gravity = @(p) 1./((0.50 - p).^2 + 0.25^2).^(3/2);
gravity2 = @(p) 1./((0.50 - p).^2 + 0.25^2).^(4/2);
gravity3 = @(p) 1./((0.50 - p).^2 + 0.50^2).^(3/2);

 gravity = @(p)  gravity(p)/integral(gravity ,-1000,1000);
gravity2 = @(p) gravity2(p)/integral(gravity2,-1000,1000);
gravity3 = @(p) gravity3(p)/integral(gravity3,-1000,1000);

figure(1)
hold on
    xlim([0 1]);
    fplot(gravity,'LineWidth',2,'Color','blue')
    fplot(gravity2,'LineWidth',2,'Color','red')
    fplot(gravity3,'LineWidth',2,'Color','magenta')
    legend('$\alpha = 3$, $\beta = 0.25$','$\alpha = 4$, $\beta = 0.25$','$\alpha = 3$, $\beta = 0.50$', 'Interpreter','LaTeX','FontSize',12)

    set(gcf,'position',[700,250,600,300])
    set(gcf,'PaperOrientation','landscape');
    exportgraphics(gcf,'../figures/gravityplot.png')
hold off
