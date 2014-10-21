%% FIXED PARAMS
N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
alpha  = .8;
t      = 1:20;

params = struct('N', N, 'alpha', alpha, 't', t);
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

ALL_MOVES  = {'pairwise','shuffle','mergesplit'};
moves      = {1,[1,2],[1,3],[1,2,3],[2,3]};

nRepeats = 10;

nExperiments = numel(moves);

mcmc_nIter = 5000;

%% EXPERIMENT

LL = cell(1,nExperiments);
runtime = zeros(1,nExperiments);

for iRepeat = 1:nRepeats    
    %% DATA GENERATION
    % note, this is independent of move choices, so held fixed for all experiments
    %generate the true data
    [y, q, n, p] = sampleState(N, mu, sigma, lambda, alpha, t);
    
    %generate the naive start state
    [q_0, n_0] = naiveStateExplanation(y,N);
    
    %% MCMC
    for iExperiment = 1:nExperiments
        %run MCMC for a bit
        tic
        [~,~,LLResult] = mcmc(p,q_0,n_0,y,alpha,'iterations',mcmc_nIter,'moves',ALL_MOVES(moves{iExperiment}));
        if isempty(LL{iExperiment})
            LL{iExperiment} = LLResult;
        else
            LL{iExperiment} = LL{iExperiment} + LLResult;
        end
        runtimeResult = toc;
        runtime(iExperiment) = runtime(iExperiment) + runtimeResult;
        
        format
        fprintf('Experiment %i, repeat %i: moves {%s}, mcmc time is %.4f sec\n', iExperiment, iRepeat, strjoin(ALL_MOVES(moves{iExperiment}),','), runtimeResult);
    end
end

%% RESULTS

%convert LL to cumulative mean LL
LL = cellfun(@(ll) ll / nRepeats, LL, 'UniformOutput', false);
MLL = cellfun(@(ll) cumsum(ll) ./ (1:length(ll)), LL, 'UniformOutput', false);

%note: display routines may be hardcoded to nExperiments...
figure(1)
hold all

legendNames = cell(1,nExperiments);
for iExperiment = 1:nExperiments
    plot(-MLL{iExperiment})
    legendNames{iExperiment} = strjoin(ALL_MOVES(moves{iExperiment}),',');
end
legend(legendNames{:})
xlabel('Iteration')
ylabel('Cumulative Mean NLL')
title('Effect of move pool on burn in time')