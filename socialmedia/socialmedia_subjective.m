clear all
close all

xpts = linspace(-1,1,10000);
bw = 0.03;

T = 100;
N = 100;
M = 3;

alpha = 2.7;
s = 0.24;
lambda = .00;
m = 0.20;

tribes = zeros(3,T);

rng(2)
sigma = 0.025;
epsilon = sigma*randn(N,T);
rng(2+1)
zeta = 0.15;
nu = zeta*randn(N,T);

rng(2)
agents(:,1) = unifrnd(-1,1,N,1);
magents = agents(:,1);


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

        %gravity equation
        gravweights(:,n,t) = (1 ./ ((agents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
        
        %properly scaled weights
        weights(:,n,t) = (gravweights(:,n,t) / (sum(gravweights(:,n,t))));
        f(n,t) = agents(n,t-1) + nu(n,t);
        %new weighted position
        agents(n,t+1) = (1-lambda)*((1-m)*sum(agents(:,t).*weights(:,n,t)) + m*f(n,t)) ...
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
    
[~,ttc ] = max(diff(consensus));
disp(['Time to Consensus: ' num2str(ttc)])


%%
close(figure(1))
figure(1)
hold on
%   ylim([-1,1])
    xlim([1,T])
    title('Social Media, Subjective (Single Sim)','FontSize',10)
    
    plot(agents')
    %plot(magents','LineWidth',6,'Color','black')
    %plot(magents','LineWidth',3,'Color','magenta')
    
    set(gcf,'position',[700,250,400,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/socialmedia_subjective.png')
hold off

