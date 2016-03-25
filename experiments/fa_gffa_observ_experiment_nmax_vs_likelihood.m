function ll = fa_gffa_observ_experiment_nmax_vs_likelihood(varargin)
	DEFAULT_NMAX   = 10:10:300;
	DEFAULT_K      = 3;
	DEFAULT_GAMMA0 = 100;
	DEFAULT_ALPHA0 = 0.35;
	DEFAULT_SILENT = false;

	parser = inputParser;
	addOptional  (parser, 'n_max',   DEFAULT_NMAX)
	addParamValue(parser, 'K',       DEFAULT_K)
	addParamValue(parser, 'gamma_0', DEFAULT_GAMMA0)
	addParamValue(parser, 'alpha_0', DEFAULT_ALPHA0)
	addParamValue(parser, 'silent',  DEFAULT_SILENT)

	parse(parser, varargin{:})
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

	%sample observations
	n_true = poissrnd(intervalGamma);
	y      = binornd(n_true, observAlpha);

	%get likelihood w/gffa
	ll_gffa = gf_forward(y, intervalGamma, alpha_0, intervalDelta);

	%get true likelihood, fa ll with each setting of n_max
	ll_fa = zeros(nExperiments,1);
	ll_tr = zeros(nExperiments,1);
	for iExperiment = 1:nExperiments
		n_max_experiment = n_max(iExperiment);

		fprintf('Experiment %d/%d\n', iExperiment, nExperiments)

		%true likelihood
		ll = zeros(n_max_experiment+1,1);
		for n = 0:n_max_experiment
			ll(n+1) = poisspdf(n, gamma_0) * prod(binopdf(y, n, alpha_0));
		end
		ll_tr(iExperiment) = sum(ll);

		%fa likelihood
		[~,~,ll] = forward_messages(y, intervalGamma, intervalDelta, observAlpha, n_max_experiment);
		ll_fa(iExperiment) = exp(ll);
	end

	if ~silent
		fprintf('nMax\ttrue\tFA\tGFFA\n')
		for iExperiment = 1:nExperiments
			fprintf('%d\t%.2f\t%.2f\t%.2f\n', n_max(iExperiment), ll_tr(iExperiment), ll_fa(iExperiment), ll_gffa)
		end
	end

	if false
		figure
		hold on
		plot(n_max, runtimeFA, 'o-', n_max, runtimeGFFA .* ones(nExperiments, 1), '-')
		title('n_{max} vs runtime, observation test')
		xlabel('n_{max}')
		ylabel('mean runtime (s)')
	end

	ll = [ll_tr, ll_fa, ll_gffa .* ones(nExperiments, 1)];
end

