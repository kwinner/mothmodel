%% FIXED PARAMS
N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
alpha  = 0.25;
t      = 1:20;

params = struct('N', N, 'alpha', alpha, 't', t);
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

ALL_MOVES  = {'pair','shuffle','cycle','mergesplit'};
% moves      = {1,[1,2],[1,4],[1,2,4],[2,4]};
% moves = {1,3,[1,3],[1,2,3],[1,2,3,4]};
% moves = {1,[1,2],3,[1,3],[1,2,3,4]};
moves = {1};

nRepeats = 1;

nExperiments = numel(moves);

mcmc_nIter = 10000;

%% EXPERIMENT

LL = cell(1,nExperiments);
runtime = zeros(1,nExperiments);

for iRepeat = 1:nRepeats    
    %% DATA GENERATION
    % note, this is independent of move choices, so held fixed for all experiments
    %generate the true data
    [y, state] = sampleState(theta, params);
    
    %generate the naive start state
    [state_0] = naiveStateExplanation(y,N);
    state_0.p = ppdf(theta, params);
    
    %% MCMC
    for iExperiment = 1:nExperiments
        %run MCMC for a bit
        tic
        state_history = mcmc(y,state_0,params,'iterations',mcmc_nIter,'moves',ALL_MOVES(moves{iExperiment}));
        if isempty(LL{iExperiment})
            LL{iExperiment} = arrayfun(@(x) loglikelihood(x,y,params), state_history);
        else
            LL{iExperiment} = LL{iExperiment} + arrayfun(@(x) loglikelihood(x,y,params), state_history);
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
figure
hold all

legendNames = cell(1,nExperiments);
for iExperiment = 1:nExperiments
    plot(-MLL{iExperiment},'linewidth',3)
    legendNames{iExperiment} = strjoin(ALL_MOVES(moves{iExperiment}),',');
end
legend(legendNames{:})
xlabel('Iteration')
ylabel('Cumulative Mean NLL')
title(sprintf('Effect of move pool on convergence, alpha = %.2f',alpha))