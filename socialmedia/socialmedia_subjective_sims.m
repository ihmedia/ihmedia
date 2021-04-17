clear all
close all

xpts = linspace(-1,1,10000);
bw = 0.05;
mph = 0.25;
mpd = 800;
sims = 1000;
extra = 100;

T = 100;
N = 100;
M = 3;

alpha = 2.7;
s = .24;
lambda = 0;
m = 0.20;

sigma = 0.025;
zeta = 0.15;

wb = parwaitbar(sims,'BarLength',10);
parfor j = 1 : sims  + extra
    j;
    agents = [];
    magents = [];
    gravweights = [];
    mgravweights = [];
    weights = [];
    mweights = [];
    kdvector = [];
    a = [];
    b = [];
    tribes = [];
    numtribes = zeros(T,1);
    f = [];   
    
    rng(j)
    epsilon = sigma*randn(N,T);
    nu = zeta*randn(N,T);
       
    agents(:,1) = unifrnd(-1,1,N,1);
    magents(:,[1:T]) = NaN(M,T);
    
	broken = 0;

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
        [b, d] = findpeaks(a(:,t)','MinPeakHeight',mph,'MinPeakDistance',mpd); 
        d = d + (d == 5000);
        numtribes(t) = length(b);
        
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

    for t = tbar + 1 : T
        if broken == 1
            magents = NaN;
            agents = NaN;
            continue
        end
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
        [b, d] = findpeaks(a(:,t)','MinPeakHeight',mph,'MinPeakDistance',mpd);
        d = d + (d == 5000);
        numtribes(t) = length(b);
                 
       
        if numtribes(t) > 1 && numtribes(t-1) == 1
            numtribes = [];   
            disp(['Simulation ' num2str(j),' is broken!'])
            broken = 1;
            continue
        end      

        if length(b) == 3
            tribes(:,t) = ((d - length(xpts)/2)/(length(xpts)/2))';
        elseif length(b) == 2
            tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0];
        elseif length(b) == 1
            tribes(:,t) = [((d - length(xpts)/2)/(length(xpts)/2))'; 0; 0];
            
        end 
  
    end
    
    if broken ~= 1
        tribecount(:,j) = numtribes;
        works(j) = 1;
    else
        works(j) = 0;
    end
    
    
    pause(rand);
    wb.progress();
end

tribecount( :, all(~tribecount,1) ) = [];

% number that reach convergence
disp(['Sims to Reach Consensus Before t = 100: ', num2str(sum(tribecount(100,[1:sims]) == 1))])


%%
close(figure(1))
figure(1)
hold on    
    ylim([0 10])
    title('Social Media, Subjective (1000 Sims)','FontSize',10)

    plot(tribecount(:,[1:sims]))     
    
    set(gcf,'position',[700,250,400,300])
    set(gcf,'PaperOrientation','landscape');
    %exportgraphics(gcf,'../figures/socialmedia_subjective_sims.png')
hold off
