function runtime = fa_gffa_observ_experiment_nmax_vs_runtime(varargin)
	DEFAULT_NITER  = 5;
	DEFAULT_NMAX   = 10:5:200;
	DEFAULT_K      = 10;
	DEFAULT_GAMMA0 = 8;
	DEFAULT_ALPHA0 = 0.5;
	DEFAULT_SILENT = false;

	parser = inputParser;
	addOptional  (parser, 'nIter',  DEFAULT_NITER)
	addOptional  (parser, 'n_max',   DEFAULT_NMAX)
	addParamValue(parser, 'K',       DEFAULT_K)
	addParamValue(parser, 'gamma_0', DEFAULT_GAMMA0)
	addParamValue(parser, 'alpha_0', DEFAULT_ALPHA0)
	addParamValue(parser, 'silent',  DEFAULT_SILENT)

	parse(parser, varargin{:})
	nIter   = parser.Results.nIter;
	n_max   = parser.Results.n_max;
	K       = parser.Results.K;
	gamma_0 = parser.Results.gamma_0;
	alpha_0 = parser.Results.alpha_0;
	silent  = parser.Results.silent;

	nExperiments = numel(n_max);

	%populate per-interval params
	intervalGamma = [gamma_0, zeros(1,K-1)];
	intervalDelta = ones(1,K-1);
	observAlpha   = alpha_0 .* ones(1,K);

	runtimeFA   = zeros(nExperiments, 1);
	runtimeGFFA = 0;
	for iter = 1:nIter
		if ~silent, fprintf('Iteration %d of %d', iter, nIter); end

		%sample observations
		n_true = poissrnd(intervalGamma);
		y      = binornd(n_true, observAlpha);
		
		%test gffa
		tStart = tic;
		ll_gffa = gf_forward(y, intervalGamma, alpha_0, intervalDelta);
		runtimeGFFA = runtimeGFFA + toc(tStart);

		%test fa with each setting of n_max
		for iExperiment = 1:nExperiments
			tStart = tic;
			[~,~,ll_fa] = forward_messages(y, intervalGamma, intervalDelta, observAlpha, n_max(iExperiment));
			% ll_fa = exp(ll_fa);
			runtimeFA(iExperiment) = runtimeFA(iExperiment) + toc(tStart);

		end

		if ~silent, fprintf(': %.2fm (~%.2fm remaining)\n', sum([runtimeFA; runtimeGFFA])/60, sum([runtimeFA; runtimeGFFA] ./ iter) * (nIter - iter)/60); end
	end

	runtimeFA   = runtimeFA   ./ nIter;
	runtimeGFFA = runtimeGFFA /  nIter;

	if ~silent
		fprintf('nMax\tFA\tGFFA\n')
		for iExperiment = 1:nExperiments
			fprintf('%d\t%.2fs\t%.2fs\n', n_max(iExperiment), runtimeFA(iExperiment), runtimeGFFA)
		end
	end

	if ~silent
		figure
		hold on
		plot(n_max, runtimeFA, 'o-', n_max, runtimeGFFA .* ones(nExperiments, 1), '-')
		title('n_{max} vs runtime, observation test')
		xlabel('n_{max}')
		ylabel('mean runtime (s)')
	end

	runtime = [runtimeFA, runtimeGFFA .* ones(nExperiments, 1)];
end

