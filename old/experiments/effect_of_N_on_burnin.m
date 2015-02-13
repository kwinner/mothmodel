%% FIXED PARAMS
%N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
alpha  = .5;
t      = 1:20;

nExperiments = 7;

mcmc_nIter = 250;

%% EXPERIMENT
LL      = cell(1,nExperiments);
runtime = zeros(1,nExperiments);
for iExperiment = 1:nExperiments
    N(iExperiment) = round(10^((iExperiment+1)/2));
    
    %generate the true data
    [y, q, n, p] = sampleState(N(iExperiment), mu, sigma, lambda, alpha, t);
    
    %generate the naive start state
    [q_0, n_0] = naiveStateExplanation(y,N(iExperiment));
    
    %run MCMC for a bit
    tic
    [~,~,LL{iExperiment}] = mcmc(p,q_0,n_0,y,alpha,'iterations',mcmc_nIter);
    runtime(iExperiment) = toc;
    
    format
    fprintf('Experiment %i: %i individuals, mcmc time is %.4f sec\n', iExperiment, N(iExperiment), runtime(iExperiment));
end

%% RESULTS

%convert LL to cumulative mean LL
MLL = cellfun(@(ll) cumsum(ll) ./ (1:length(ll)), LL, 'UniformOutput', false);

%note: display routines may be hardcoded to nExperiments...
figure(1)
plot(N, runtime,'o-')
title('#indiv vs runtime')
xlabel('N')
ylabel('runtime')

figure(2)
hold on

for iExperiment = 1:nExperiments
    subplot(nExperiments, 1, iExperiment)
    plot(-MLL{iExperiment})
    title(num2str(N(iExperiment)))
end