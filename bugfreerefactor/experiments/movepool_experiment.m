function NLL_avg_cumavg = movepool_experiment( mu, sigma, lambda, N, alpha, T, varargin )

DEFAULT_MCMC_ITERATIONS = 2500;
DEFAULT_REPEATS         = 10;
DEFAULT_USE_ARS         = true;

parser = inputParser;
addParamValue(parser, 'nMCMCIterations', DEFAULT_MCMC_ITERATIONS);
addParamValue(parser, 'nRepeats',        DEFAULT_REPEATS);
addParamValue(parser, 'useARS',          DEFAULT_USE_ARS);

parser.parse(varargin{:});
nMCMCIterations = parser.Results.nMCMCIterations;
nRepeats        = parser.Results.nRepeats;
useARS          = parser.Results.useARS;

%setup movepools: the sets of moves to experiment over
ALL_MOVES = {'pair', 'shuffle', 'cycle', 'mergesplit'};
movepools = {1, 3, [1,3]};
colors = 'rgb';

nExperiments = numel(movepools);

%LL = record for LL of all experiments
LL = zeros(nRepeats, nExperiments, nMCMCIterations);
for iRepeat = 1:nRepeats
	%generate the true observations
	[y, n, S, Z, Q] = samplestate(mu, sigma, lambda, N, alpha, T);
	P = birthdeath_pmf(mu, sigma, lambda, N, T);

	%generate the canonical start state
	Q_0 = canonical_state(y,N);

	%run MCMC
	for iExperiment = 1:nExperiments
		tic
		[ ~, Q_history] = mcmc(y, Q_0, P, N, alpha, T, 'nIterations', nMCMCIterations, 'moves', ALL_MOVES(movepools{iExperiment}));

		abundance_history = cellfun(@abundance, Q_history, 'UniformOutput', false);
		LL(iRepeat, iExperiment,:) = cellfun(@(q,n) latentvar_LL(q, P) + latentvar_obs_LL(y, n, alpha), Q_history, abundance_history);
		runtime = toc;
		fprintf('Experiment %i, repeat %i: moves {%s}, time is %.4fs\n', iExperiment, iRepeat, strjoin(ALL_MOVES(movepools{iExperiment}),','), runtime);
	end
end

LL_avg = squeeze(mean(LL,1));
NLL_avg_cumavg = -1*cumsum(LL_avg,2)./repmat(1:nMCMCIterations,nExperiments,1);

figure; hold on
legendNames = cell(1,nExperiments);
for iExperiment = 1:nExperiments
	plot(NLL_avg_cumavg(iExperiment,:), colors(iExperiment));
	legendNames{iExperiment} = strjoin(ALL_MOVES(movepools{iExperiment}), ',');
end
legend(legendNames{:})
xlabel('Iteration')
ylabel('Cumulative Mean NLL')
title(['Effect of move pool on convergence, alpha = ', num2str(alpha)])

keyboard

end

