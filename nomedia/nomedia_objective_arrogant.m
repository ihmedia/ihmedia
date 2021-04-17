clear all
close all

xpts = linspace(-1,1,10000);
bw = 0.05;
mph = 0.25;
mpd = 800;

T = 100;
N = 100;
M = 3;

alpha = 8;
s = .25;
m = 0;
lambda = .01917;

rng(987)
sigma = 0.025;
epsilon = sigma*randn(N,T);
 
tribes = zeros(3,T);

%agents(:,1) = trnd(7,N,1);
rng(987) 
agents(:,1) = unifrnd(-1,1,N,1);
magents(:,[1:T]) = NaN(M,T);


%%
for t = 1 : T
    for n = 1 : N     
        %gravity equation
        gravweights(:,n,t) = (1 ./ ((agents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);

        %properly scaled weights
        weights(:,n,t) = gravweights(:,n,t) / sum(gravweights(:,n,t));

        %new weighted position
        agents(n,t+1) = (1-lambda)*sum(agents(:,t).*weights(:,n,t)) + epsilon(n,t);
    end

    [a(:,t),~] = ksdensity(agents(:,t),xpts,'Bandwidth',bw);
    kdvector = a(:,t)';
    [b, d] = findpeaks(a(:,t)','MinPeakHeight',0.20,'MinPeakHeight',0.20,'MinPeakDistance',800);
     
    if length(b) == 2
        tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0];
    elseif length(b) == 1
        tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0; 0];
    elseif length(b) == 3
        tribes(:,t) = ((d - length(xpts)/2)/(length(xpts)/2))';
    end
    
    %if only one tribe, consensus reached
    if length(b) == 1
        consensus(t) = 1;
    end
    
    if length(b) <= 3
        break
    end
end
    
tbar = t;

%%    
for t = tbar + 1 : T
    for n = 1 : N     
        magents(:,t) = tribes(:,t-1);

        %gravity equation
        gravweights(:,n,t) = (1 ./ ((agents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
        mgravweights(:,n,t) = (1 ./ ((magents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);

        %properly scaled weights
        weights(:,n,t) = (1-m)*(gravweights(:,n,t) / (sum(gravweights(:,n,t))));
        mweights(:,n,t) = (mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));
        mweights(:,n,t) = m*(mweights(:,n,t).*(magents(:,t) ~= 0)/sum(mweights(:,n,t).*(magents(:,t) ~= 0)));
        %mweights(:,n,t) = mweights(:,n,t).*(magents(:,t) ~= 0) / sum(mweights(:,n,t).*(magents(:,t) ~= 0));

        %new weighted position
        agents(n,t+1) = (1-lambda)*( sum(agents(:,t).*weights(:,n,t))  ...
                        + sum(magents(:,t).*mweights(:,n,t))) ...
                        + epsilon(n,t);
    end
      
    %calculate tribes
    [a(:,t),~] = ksdensity(agents(:,t),xpts,'Bandwidth',bw);
    kdvector = a(:,t)';
    [b, d] = findpeaks(a(:,t)','MinPeakHeight',0.20,'MinPeakDistance',800);
    
    numtribes(t) = length(b);
    
    if length(b) == 3
        tribes(:,t) = ((d - length(xpts)/2)/(length(xpts)/2))';
    elseif length(b) == 2
        tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0];
    elseif length(b) == 1
        tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0; 0];
    end
	        
    %if only one tribe, consensus reached
    if length(b) == 1
        consensus(t) = 1;
    end  
end
    
magents(:,t) = tribes(:,t-1);


%%
close(figure(1))
figure(1)
hold on
%   ylim([-1,1])
    xlim([1,T])
    
    title('Objective Dynamics, Low IH (Single Sim)','FontSize',10)
    plot(agents')
    %plot(magents','LineWidth',6,'Color','black')
    %plot(magents','LineWidth',3,'Color','magenta')
    set(gcf,'position',[700,250,400,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/nomedia_objective_arrogant.png')
hold off


%%
close(figure(2))
figure(2)
hold on
    set(gca,'ytick',[]);
    set(gca,'DefaultLineLineWidth',2)
    xticks([1,2500,5000,7500,10000])
    xticklabels({'-1','-0.5','0','0.5','1'})
    title('Four Tribes: t = 30','FontSize',10)
    
    findpeaks(a(:,20)','MinPeakHeight',mph,'MinPeakDistance',mpd)
    set(gcf,'position',[700,250,315,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/nomedia_objective_arrogant_kdensity1.png')
hold off


%%
close(figure(3))
figure(3)
hold on
    set(gca,'ytick',[]);
    set(gca,'DefaultLineLineWidth',2)
    xticks([1,2500,5000,7500,10000])
    xticklabels({'-1','-0.5','0','0.5','1'})
    title('Two Tribes: t = 60','FontSize',10)
    
    findpeaks(a(:,60)','MinPeakHeight',mph,'MinPeakDistance',mpd)
    set(gcf,'position',[700,250,315,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/nomedia_objective_arrogant_kdensity2.png')
hold off


%%
close(figure(4))
figure(4)
hold on
    set(gca,'ytick',[]);
    set(gca,'DefaultLineLineWidth',2)
    xticks([1,2500,5000,7500,10000])
    xticklabels({'-1','-0.5','0','0.5','1'})
    title('Consensus: t = 90','FontSize',10)
    
    findpeaks(a(:,80)','MinPeakHeight',mph,'MinPeakDistance',mpd)
    
    set(gcf,'position',[700,250,315,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/nomedia_objective_arrogant_kdensity3.png')
hold off
