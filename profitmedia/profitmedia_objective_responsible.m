clear all
close all

xpts = linspace(-1,1,10000);
bw = 0.03;

T = 100;
N = 100;
M = 3;

alpha = 4;
s = .25;
lambda = .01917;
m = .49;

rng(987)
sigma = 0.025;
epsilon = sigma*randn(N,T);

rng(987+1)
zeta = 0.005;
nu = zeta*randn(T,1);

tribes = zeros(3,T);

%agents(:,1) = trnd(7,N,1);
rng(987)
agents(:,1) = unifrnd(-1,1,N,1);
magents(:,[1:T]) = NaN(M,T);

true_consensus = zeros(T,1);


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
    
    if length(b) <= 3
        break
    end
end
    
tbar = t;

%%    
for t = tbar + 1 : T
    for n = 1 : N     
        magents([1:3],t) = tribes(:,t-1);
        magents(4,t) = 0 + nu(t);      

        %gravity equation
        gravweights(:,n,t) = (1 ./ ((agents(:,t) - agents(n,t)).^2 + s^2 ).^alpha/2);

        %properly scaled weights
        weights(:,n,t) = (1-m)*(gravweights(:,n,t) / (sum(gravweights(:,n,t))));
        
        numprofmedia = nnz(tribes([1:3],t-1));
        
        if numprofmedia == 3
            mgravweights([1:3],n,t) = (1 ./ ((magents([1:3],t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mgravweights(4,n,t) = (1 ./ ((magents(4,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mweights(:,n,t) = m*(mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));
        elseif numprofmedia == 2
            mgravweights([1:2],n,t) = (1 ./ ((magents([1:2],t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mgravweights(3,n,t) = 0;    
            mgravweights(4,n,t) = (1 ./ ((magents(4,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mweights(:,n,t) = m*(mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));    
        elseif numprofmedia == 1 && abs(mean(magents(1,[t-2:t]))) > 0.025
            mgravweights(1,n,t) = (1 ./ ((magents(1,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mgravweights([2:3],n,t) = 0;    
            mgravweights(4,n,t) = (1 ./ ((magents(4,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mweights(:,n,t) = m*(mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));    
        else
            %mgravweights(1,n,t) = (1 ./ ((magents(1,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            magents(1,t) = 0;
            mgravweights([1:3],n,t) = 0;    
            mgravweights(4,n,t) = (1 ./ ((magents(4,t) - agents(n,t)).^2 + s^2 ).^alpha/2);
            mweights(:,n,t) = m*(mgravweights(:,n,t) / (sum(mgravweights(:,n,t))));                
        end

        %new weighted position
        agents(n,t+1) = (1-lambda)*( sum(agents(:,t).*weights(:,n,t))  ...
                        + sum( magents(:,t).*mweights(:,n,t)) ) ...
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
    if length(b) == 1 && abs(mean(magents(1,[t-2:t]))) < 0.025
        true_consensus(t) = 1;
        abs(   mean(magents(1,[tbar+1:t]))      );
    end  
end
    


%%
merger = (diff(numtribes) == 0);
merger([T:T+2]) = 1;

magents(magents==0) = NaN;

[y,z] = max(diff(magents(2,:)));
magents(2,z) = NaN;
magents(2,z-1) = magents(2,z);
magents(3,z+1) = magents(2,z+1);


%%
close(figure(1))
figure(1)
hold on
%   ylim([-1,1])
    xlim([1,T])
    title('Profit & Responsible Media, Objective','FontSize',10)
    
    plot(agents')
    plot(magents([1:3],:)','LineWidth',6,'Color','black')
    plot(magents([1:3],:)','LineWidth',3,'Color','green')
    plot(magents([4],:)','LineWidth',6,'Color','black')
    plot(magents([4],:)','LineWidth',3,'Color','magenta')
    
    set(gcf,'position',[700,250,400,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/profitmedia_objective_responsible.png')
hold off

