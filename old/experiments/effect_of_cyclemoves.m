% FIXED PARAMS
N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
alpha  = 1;
t      = 1:20;

% gross
% N      = 73;
% mu     = 9.5;
% sigma  = 2.8;
% lambda = 1.488;
% alpha  = .27;
% t      = 1:20;


params = struct('N', N, 'alpha', alpha, 't', t);
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

ALL_MOVES  = {'pair','shuffle','cycle', 'mergesplit'};
moves      = {1,3,[1,3],[2,3],[1,4],[1,2,4]};
colors = 'bgrcmy';

nRepeats = 10;

nExperiments = numel(moves);

mcmc_nIter = 2000;

use_ARS = false;

%% EXPERIMENT

LL = zeros(nRepeats,nExperiments,mcmc_nIter);
runtime = zeros(nRepeats,nExperiments);

for iRepeat = 1:nRepeats    
    %% DATA GENERATION
    % note, this is independent of move choices, so held fixed for all experiments
    
    %generate the observations
    y = sampleState(theta, params);
    
    %generate the naive start state
    state_0 = naiveStateExplanation(y,N);
    state_0.p = ppdf(theta, params);
    
    %% MCMC
    for iExperiment = 1:nExperiments
        %run MCMC for a bit
        tic
        state = mcmc(y,state_0,params,'iterations',mcmc_nIter,'moves',ALL_MOVES(moves{iExperiment}),'use_ARS',use_ARS);

        LL(iRepeat,iExperiment,:) = LL(iRepeat,iExperiment,:) + shiftdim(arrayfun(@(s) loglikelihood(s, y, params), state),-1);
     
        runtime(iRepeat,iExperiment) = toc;     
        
        format
        fprintf('Experiment %i, repeat %i: moves {%s}, mcmc time is %.4f sec\n', iExperiment, iRepeat, strjoin(ALL_MOVES(moves{iExperiment}),','), runtime(iRepeat,iExperiment));
    end
end

% calculate average and cumavg
runtime_avg = squeeze(mean(runtime,1));
LL_avg = squeeze(mean(LL,1));
NLL_avg_cumavg = -1*cumsum(LL_avg,2)./repmat(1:mcmc_nIter,nExperiments,1);

% plot
figure; hold on
legendNames = cell(1,nExperiments);
for iExperiment = 1:nExperiments
    plot(NLL_avg_cumavg(iExperiment,:), colors(iExperiment), 'linewidth', 3);
    legendNames{iExperiment} = strjoin(ALL_MOVES(moves{iExperiment}),',');
end
legend(legendNames{:})
xlabel('Iteration')
ylabel('Cumulative Mean NLL')
title(sprintf('Effect of move pool on convergence, alpha = %d',alpha))
