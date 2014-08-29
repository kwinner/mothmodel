N_EXPERIMENTS = 6;
COLORS        = {'r','g','b','c','m','y'};
figure
hold on
for iExperiment = 1:N_EXPERIMENTS
    N_ITER = 1000;
    t      = 1:10;
    T      = numel(t);
    NQ     = (T+1)*(T+2)/2;
    N = 10^iExperiment;
%     N      = 100;
    mu     = 5;
    sigma  = 2.5;
    lambda = 2;
%     alpha  = (iExperiment-1)/(N_EXPERIMENTS-1);
%     if iExperiment == N_EXPERIMENTS
%         alpha = .9999999;
%     end
alpha = .5;
    epsilon = .01;

    q = cell(1,N_ITER);
    n = cell(1,N_ITER);
    normalizedLikelihood1 = cell(1,N_ITER);
    
    [y, ~, ~, p] = generateMothData(N, mu, sigma, lambda, alpha, t);
    
    %old naive generation method
%     n{1} = y;
% %     n{1}(1) = N-sum(y);
%     q{1} = zeros(T+1);
%     for i = 1:T
%         q{1}(i,i+1) = n{1}(i);
%     end
%     q{1}(1,1) = q{1}(1,1) + N - sum(sum(q{1}));

    %new, individual, *correct* naive generation method
    nAlive = 0;
    nDead  = 0;
    births = zeros(1,T+1);
    deaths = zeros(1,T+1);
    q{1}   = zeros(T+1,T+1);
    for i = 1:T
        if y(i) >= nAlive %more individuals needed to explain y(i)
            %birth(I(i-1),y(i)-#alive)
            babies = y(i) - nAlive;
            births(i) = births(i) + babies;
            nAlive = nAlive + babies;
        else %too many individuals alive, kill some
            %kill(I(i-1),y(i)-nAlive)
            quota = nAlive - y(i);
            for j = 1:i
                if deaths(j) < births(j)
                    avail = births(j) - deaths(j); %how many we can kill from j
                    victims = min(quota,avail);        %how many we actually need to kill
                    %kill them:
                    deaths(j) = deaths(j) + victims;
                    nAlive    = nAlive - victims;
                    nDead     = nDead + victims;
                    quota     = quota - victims;
                    q{1}(j,i) = q{1}(j,i) + victims;
                    %do we need to keep killing??
                    if quota <= 0; break; end
                end
            end
        end
    end
    nObserved = nAlive + nDead;
    %kill the survivors
    for i = 1:T+1
        if deaths(i) < births(i)
            victims = births(i) - deaths(i);
            q{1}(i,T+1) = q{1}(i,T+1) + victims;
            nAlive = nAlive - victims;
            nDead  = nDead + victims;
            if nAlive <= 0; break; end
        end
    end
    %throw the 'unobserved' individuals in at the beginning: q(1,1)
    q{1}(1,1) = N - sum(q{1}(:));
    n{1} = abundancy(q{1});
    
    i = 1;
    LL = zeros(N_ITER,1);
    LL(i,:) = mothLikelihood(p,q{1},n{1},y,alpha);
    
    i = 2;
    while i<=N_ITER
        DELTA = randsample(NQ,2,false);
        [a,b] = ind2sub_triu(T+1, DELTA(1));
        [c,d] = ind2sub_triu(T+1, DELTA(2));
        
        [q{i}, n{i}, normalizedLikelihood1{i}] = gibbs(p,q{i-1},n{i-1},y,alpha,a,b,c,d,'logspace',true);
        
        if q{i} == q{i-1}
            continue
        end
        
        LL(i) = mothLikelihood(p,q{i},n{i},y,alpha);
            
%         convergence = std(LL(1:i))
%         if convergence < epsilon
%             break
%         end
        
        i = i+1;
    end
    
    plot(-LL, [COLORS{iExperiment}, '-'])
    legendLabel{iExperiment} = num2str(alpha);
end

legend(legendLabel);