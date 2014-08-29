%% FIXED PARAMS
N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
alpha  = .5;
t      = 1:20;

ALL_MOVES  = {'pairwise','shuffle','mergesplit'};
moves      = {1,[1,2],[1,3],[1,2,3],[2,3]};    

nExperiments = numel(moves);

mcmc_nIter = 2500;

%% DATA GENERATION 
% note, this is independent of move choices, so held fixed for all experiments
%generate the true data
[y, q, n, p] = sampleState(N, mu, sigma, lambda, alpha, t);

%generate the naive start state
[q_0, n_0] = naiveStateExplanation(y,N);

%% EXPERIMENT
LL      = cell(1,nExperiments);
runtime = zeros(1,nExperiments);
for iExperiment = 1:nExperiments
    %run MCMC for a bit
    tic
    [~,~,LL{iExperiment}] = mcmc(p,q_0,n_0,y,alpha,'iterations',mcmc_nIter,'moves',ALL_MOVES(moves{iExperiment}));
    runtime(iExperiment) = toc;
    
    format
    fprintf('Experiment %i: moves {%s}, mcmc time is %.4f sec\n', iExperiment, strjoin(ALL_MOVES(moves{iExperiment}),','), runtime(iExperiment));
end

%% RESULTS

%note: display routines may be hardcoded to nExperiments...
figure(1)
hold all

legendNames = cell(1,nExperiments);
for iExperiment = 1:nExperiments
    plot(-LL{iExperiment})
    legendNames{iExperiment} = strjoin(ALL_MOVES(moves{iExperiment}),',');
end
legend(legendNames{:})
xlabel('Iteration')
ylabel('NLL')
title('Effect of move pool on burn in time')