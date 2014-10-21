%% FIXED PARAMS
N      = 100;
mu     = 8;
sigma  = 4;
lambda = 3;
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);
% alpha  = .5;
t      = 1:20;
moves  = {'pairwise'};
% moves  = {'pairwise','shuffle','mergesplit'};

alpha = [0.5 1];
nExperiments = numel(alpha);

mcmc_nIter = 1500;


%% EXPERIMENT
LL      = cell(1,nExperiments);
runtime = zeros(1,nExperiments);
Q = cell(1,nExperiments);
for iExperiment = 1:nExperiments
    params = struct('N', N, 'alpha', alpha(iExperiment), 't', t);
    
    %generate the observations
    y = sampleState(theta, params);
    
    %generate the naive start state
    state_0 = naiveStateExplanation(y,N);
    state_0.p = ppdf(theta, params);
    
    %run MCMC for a bit
    tic
    [Q] = mcmc(p,q_0,n_0,y,alpha(iExperiment),'iterations',mcmc_nIter,'moves',moves);
    runtime(iExperiment) = toc;
    
    format
    fprintf('Experiment %i: alpha %f, mcmc time is %.4f sec\n', iExperiment, alpha(iExperiment), runtime(iExperiment));
end

%% RESULTS

%convert LL to cumulative mean LL
MLL = cellfun(@(ll) cumsum(ll) ./ (1:length(ll)), LL, 'UniformOutput', false);

%note: display routines may be hardcoded to nExperiments...
figure(1)
plot(alpha, runtime,'o-')
title('alpha vs runtime')
xlabel('alpha')
ylabel('runtime')

figure(2)
hold on

for iExperiment = 1:nExperiments
    subplot(nExperiments, 1, iExperiment)
    plot(-MLL{iExperiment})
    title(num2str(alpha(iExperiment)))
end