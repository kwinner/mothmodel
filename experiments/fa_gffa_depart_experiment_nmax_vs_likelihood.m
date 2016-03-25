function ll = fa_gffa_depart_experiment_nmax_vs_likelihood(varargin)
	DEFAULT_NMAX    = 2:1:20;
	DEFAULT_K       = 10;
	DEFAULT_GAMMA_0 = 10;
	DEFAULT_DELTA   = 0.75;
	DEFAULT_SILENT  = false;

	parser = inputParser;
	addOptional  (parser, 'n_max',   DEFAULT_NMAX)
	addParamValue(parser, 'K',       DEFAULT_K)
	addParamValue(parser, 'gamma_0', DEFAULT_GAMMA_0)
	addParamValue(parser, 'delta',   DEFAULT_DELTA)
	addParamValue(parser, 'silent',  DEFAULT_SILENT)

	parse(parser, varargin{:})
	n_max   = parser.Results.n_max;
	K       = parser.Results.K;
	gamma_0 = parser.Results.gamma_0;
	delta   = parser.Results.delta;
	silent  = parser.Results.silent;

	nExperiments = numel(n_max);

	%populate per-interval params
	intervalGamma = [gamma_0, zeros(1,K-1)];
	intervalDelta = delta .* ones(1,K-1);
	observAlpha   = ones(1,K);

	%sample observations
	n_true = poissrnd(intervalGamma);
	for i = 2:K
		n_true(i) = binornd(n_true(i-1), intervalDelta(i-1));
	end
	y = binornd(n_true, observAlpha);

	%get likelihood w/gffa
	ll_gffa = gf_forward(y, intervalGamma, observAlpha(1), intervalDelta);

	%true likelihood
	ll_tr = poisspdf(n_true(1), gamma_0) .* prod(binopdf(n_true(2:end), n_true(1:end-1), intervalDelta));

	%get true likelihood, fa ll with each setting of n_max
	ll_fa = zeros(nExperiments,1);
	for iExperiment = 1:nExperiments
		n_max_experiment = n_max(iExperiment);

		fprintf('Experiment %d/%d\n', iExperiment, nExperiments)

		%fa likelihood
		[~,~,ll] = forward_messages(y, intervalGamma, intervalDelta, observAlpha, n_max_experiment);
		ll_fa(iExperiment) = exp(ll);
	end

	if ~silent
		fprintf('nMax\ttrue\tFA\tGFFA\n')
		for iExperiment = 1:nExperiments
			fprintf('%d\t%.2f\t%.2f\t%.2f\n', n_max(iExperiment), ll_tr, ll_fa(iExperiment), ll_gffa)
		end
	end

	if ~silent
		figure
		hold on
		plot(n_max, ll_tr .* ones(nExperiments, 1), '-', n_max(~isnan(ll_fa)), ll_fa(~isnan(ll_fa)), '-', n_max, ll_gffa .* ones(nExperiments, 1), '-')
		title('n_{max} vs likelihood, departure test')
		xlabel('n_{max}')
		ylabel('likelihood')
	end

	ll = [ll_tr .* ones(nExperiments, 1), ll_fa, ll_gffa .* ones(nExperiments, 1)];
end

