clear all
close all

xpts = linspace(-1,1,10000);
bw = 0.05;
mph = 0.25;
mpd = 800;

T = 100;
N = 100;
M = 3;

alpha = 4;
s = .25;
lambda = .01917;
m = 0.05;

bias = 0.25;

rng(987)
sigma = 0.025;
epsilon = sigma*randn(N,T);
rng(987+1)
sigmaeta = 0.025;
eta = sigmaeta*randn(2,T);

tribes = zeros(3,T);

%agents(:,1) = trnd(7,N,1);
rng(987) 
agents(:,1) = unifrnd(-1,1,N,1);
     

%%
for t = 1 : T
    
    magents(1,t) = NaN;
    magents(2,t) = NaN;        
        
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
    [b, d] = findpeaks(a(:,t)','MinPeakHeight',mph,'MinPeakDistance',mpd);   
    numtribes(t) = length(b);
        
    if length(b) <= 3
        break
    end
end
    
tbar = t;
   
for t = tbar + 1 : T  
    magents(1,t) =  bias + eta(1,t);
    magents(2,t) = -bias + eta(2,t);
        
    for n = 1 : N  

        %gravity equation
        gravweights(:,n,t) = (1 ./ ((agents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
        mgravweights(:,n,t) = (1 ./ ((magents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);

        %properly scaled weights
        weights(:,n,t) = (1-m)*(gravweights(:,n,t) / (sum(gravweights(:,n,t))));
        mweights(:,n,t) = m*(mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));
        %mweights(:,n,t) = mweights(:,n,t).*(magents(:,t) ~= 0) / sum(mweights(:,n,t).*(magents(:,t) ~= 0));

        %new weighted position
        agents(n,t+1) = (1-lambda)*( sum(agents(:,t).*weights(:,n,t))  ...
                        + sum(magents(:,t).*mweights(:,n,t))) ...
                        + epsilon(n,t);
    end
      
    %calculate tribes
    [a(:,t),~] = ksdensity(agents(:,t),xpts,'Bandwidth',bw);
    kdvector = a(:,t)';
    [b, d] = findpeaks(a(:,t)','MinPeakHeight',mph,'MinPeakDistance',mpd); 
    numtribes(t) = length(b);
   
end


%%
close(figure(1))
figure(1)
hold on
%   ylim([-1,1])
    xlim([1,T])
    title('Moderate Agenda Media (Single Sim)','FontSize',10)
    
    plot(agents')
    plot(magents(1,:)','LineWidth',6,'Color','black')
    plot(magents(1,:)','LineWidth',3,'Color','red')
    plot(magents(2,:)','LineWidth',6,'Color','black')
    plot(magents(2,:)','LineWidth',3,'Color','blue')
    
    set(gcf,'position',[700,250,400,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/agendamedia_moderate.png')
hold off

