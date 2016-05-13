function [runtimeFA, runtimeGFFA, ll_fa, ll_gffa] = fa_gffa_full_experiment_nmax(varargin)
	DEFAULT_NITER       = 5;
	DEFAULT_NMAX        = 1:1:50;
	DEFAULT_K           = 0;
	DEFAULT_MU          = 8;
	DEFAULT_SIGMA       = 4;
	DEFAULT_LAMBDA      = 3;
	DEFAULT_T           = 1:4:20;
	DEFAULT_NHAT        = 100;
	DEFAULT_GAMMA       = [];
	DEFAULT_DELTA       = [];
	DEFAULT_ALPHA       = .2;
	DEFAULT_SILENT      = false;

	parser = inputParser;
	addOptional  (parser, 'nIter',  DEFAULT_NITER)
	addOptional  (parser, 'n_max',  DEFAULT_NMAX)
	addParamValue(parser, 'K',      DEFAULT_K)
	addParamValue(parser, 'mu',     DEFAULT_MU)
	addParamValue(parser, 'sigma',  DEFAULT_SIGMA)
	addParamValue(parser, 'lambda', DEFAULT_LAMBDA)
	addParamValue(parser, 'T',      DEFAULT_T)
	addParamValue(parser, 'N_hat',  DEFAULT_NHAT)
	addParamValue(parser, 'gamma',  DEFAULT_GAMMA)
	addParamValue(parser, 'delta',  DEFAULT_DELTA)
	addParamValue(parser, 'alpha',  DEFAULT_ALPHA)
	addParamValue(parser, 'silent', DEFAULT_SILENT)

	parse(parser, varargin{:})
	nIter         = parser.Results.nIter;
	n_max         = parser.Results.n_max;
	K             = parser.Results.K;
	mu            = parser.Results.mu;
	sigma         = parser.Results.sigma;
	lambda        = parser.Results.lambda;
	T             = parser.Results.T;
	N_hat         = parser.Results.N_hat;
	intervalGamma = parser.Results.gamma;
	intervalDelta = parser.Results.delta;
	observAlpha   = parser.Results.alpha;
	silent        = parser.Results.silent;

	nExperiments = numel(n_max);

	%populate per-interval params
	if K == 0
		K = numel(T);
	end
	if numel(observAlpha) < K
		observAlpha = observAlpha(1) .* ones(1,K);
	end
	if isempty(intervalGamma == []) && isempty(intervalDelta == [])
		arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
		rateFunc     = @arrivalDistn.pdf;

		%service distribution
		serviceDistn = makedist('Exp', 'mu', lambda);

		intervalGamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
		intervalDelta = survival_prob(serviceDistn, T);
	elseif isempty(intervalGamma == []) || isempty(intervalDelta == [])
		cat('gamma and delta must both be specified or unspecified.\n')
		return
	end

	runtimeFA   = zeros(nExperiments, 1);
	runtimeGFFA = 0;
	ll_fa       = zeros(nExperiments, 1);
	ll_gffa     = 0;
	for iter = 1:nIter
		if ~silent, fprintf('Iteration %d of %d', iter, nIter); end

		%sample observations
		n_true = poissrnd(intervalGamma);
		y      = binornd(n_true, observAlpha);
		
		%test gffa
		tStart = tic;
		ll_gffa = ll_gffa + gf_forward(y, intervalGamma, observAlpha(1), intervalDelta);
		runtimeGFFA = runtimeGFFA + toc(tStart);

		%test fa with each setting of n_max
		for iExperiment = 1:nExperiments
			tStart = tic;
			[~,~,ll_temp] = forward_messages(y, intervalGamma, intervalDelta, observAlpha, n_max(iExperiment));
			% ll_temp = exp(ll_temp);
			ll_fa(iExperiment) = ll_fa(iExperiment) + ll_temp;
			runtimeFA(iExperiment) = runtimeFA(iExperiment) + toc(tStart);

		end

		if ~silent, fprintf(': %.2fm (~%.2fm remaining)\n', sum([runtimeFA; runtimeGFFA])/60, sum([runtimeFA; runtimeGFFA] ./ iter) * (nIter - iter)/60); end
	end

	runtimeFA   = runtimeFA   ./ nIter;
	runtimeGFFA = runtimeGFFA /  nIter;
	ll_fa       = ll_fa       ./ nIter;
	ll_gffa     = ll_gffa     /  nIter;

	if ~silent
		fprintf('nMax\tFA run\tPGFFA run\tFA ll\tGFFA ll\n')
		for iExperiment = 1:nExperiments
			fprintf('%d\t%.2fs\t%.2fs\t\t%.2f\t%.2f\n', n_max(iExperiment), runtimeFA(iExperiment), runtimeGFFA, ll_fa(iExperiment), ll_gffa)
		end
	end

	if ~silent
		figure
		hold on
		plot(n_max, runtimeFA, 'o-', n_max, runtimeGFFA .* ones(nExperiments, 1), '-')
		title('n_{max} vs runtime, full test')
		xlabel('n_{max}')
		ylabel('mean runtime (s)')

		figure
		hold on
		plot(n_max, ll_fa, 'o-', n_max, ll_gffa .* ones(nExperiments, 1), '-')
		title('n_{max} vs loglikelihood, full test')
		xlabel('n_{max}')
		ylabel('mean LL')
	end
end

