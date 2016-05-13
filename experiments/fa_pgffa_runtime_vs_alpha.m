%written 5/9/16
function [runtimeFA, runtimePGFFA, runtimeFAatMaxYByAlpha, runtimeFAatLLEquality, runtimeFAatTwiceNHat] = fa_pgffa_runtime_vs_alpha(varargin)
	DEFAULT_NITER		= 10;
	DEFAULT_ALPHA       = (0.05:.05:1)';

	DEFAULT_K           = 0;
	DEFAULT_MU          = 8;
	DEFAULT_SIGMA       = 4;
	DEFAULT_LAMBDA      = 3;
	DEFAULT_T           = 1:4:20;
	DEFAULT_N_HAT       = 10;
	DEFAULT_N_MAX		= 1:2:DEFAULT_N_HAT * 1.5 + 1;
	DEFAULT_GAMMA       = [];
	DEFAULT_DELTA       = [];
	DEFAULT_SILENT      = false;
	DEFAULT_LL_EPSILON  = 0.001;

	parser = inputParser;
	addParamValue(parser, 'nIter',   DEFAULT_NITER)
	addParamValue(parser, 'alpha',   DEFAULT_ALPHA)
	addParamValue(parser, 'K',       DEFAULT_K)
	addParamValue(parser, 'mu',      DEFAULT_MU)
	addParamValue(parser, 'sigma',   DEFAULT_SIGMA)
	addParamValue(parser, 'lambda',  DEFAULT_LAMBDA)
	addParamValue(parser, 'T',       DEFAULT_T)
	addParamValue(parser, 'N_hat',   DEFAULT_N_HAT)
	addParamValue(parser, 'n_max',   DEFAULT_N_MAX)
	addParamValue(parser, 'gamma',   DEFAULT_GAMMA)
	addParamValue(parser, 'delta',   DEFAULT_DELTA)
	addParamValue(parser, 'silent',  DEFAULT_SILENT)
	addParamValue(parser, 'epsilon', DEFAULT_LL_EPSILON)

	parse(parser, varargin{:})
	nIter         = parser.Results.nIter;
	observAlpha   = parser.Results.alpha;
	K             = parser.Results.K;
	mu            = parser.Results.mu;
	sigma         = parser.Results.sigma;
	lambda        = parser.Results.lambda;
	T             = parser.Results.T;
	N_hat         = parser.Results.N_hat;
	n_max         = parser.Results.n_max;
	intervalGamma = parser.Results.gamma;
	intervalDelta = parser.Results.delta;
	silent        = parser.Results.silent;
	llEpsilon     = parser.Results.epsilon;

	%populate per-interval params
	if K == 0
		K = numel(T);
	end
	%fa supports per-sample values for alpha, but we are assuming a constant detection rate
	%so alpha needs to be expanded
	if size(observAlpha, 1) == 1 && size(observAlpha, 2) > 1 && silent ~= true
		warning('Warning: alpha should probably be a column vector. dim 2 is used for per-sample alphas')
	end
	if size(observAlpha, 2) == 1
		observAlpha = repmat(observAlpha, 1, K);
	end
	if isempty(intervalGamma == []) && isempty(intervalDelta == [])
		arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
		rateFunc     = @arrivalDistn.pdf;

		%service distribution
		serviceDistn = makedist('Exp', 'mu', lambda);

		intervalGamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
		intervalDelta = survival_prob(serviceDistn, T);
	elseif isempty(intervalGamma == []) || isempty(intervalDelta == [])
		error('gamma and delta must both be specified or unspecified.')
	end
	if max(n_max) < 1.5 * N_hat
		warning('Warning: maximum value of n_hat should be at least twice N_hat, added manually.');
		n_max = [n_max, ceil(N_hat * 1.5)];
	end

	nAlpha = size(observAlpha, 1);
	nThresh = numel(n_max);

	runtimeFA              = zeros(nAlpha, nThresh);
	runtimePGFFA           = zeros(nAlpha, 1);
	runtimeFAatMaxYByAlpha = zeros(nAlpha, 1); %has to be tracked separately because corresponding iThresh depends on ymax
	runtimeFAatLLEquality  = zeros(nAlpha, 1); %also would depend on y
	llFA		           = zeros(nAlpha, nThresh);
	llPGFFA                = zeros(nAlpha, 1);
	for iAlpha = 1:nAlpha
		if ~silent, fprintf('---Experiment %d of %d, alpha = %0.2f---\n', iAlpha, nAlpha, observAlpha(iAlpha)); end
		for iter = 1:nIter
			if ~silent, fprintf('Iteration %d of %d', iter, nIter); end

			%sample observations
			n_true = poissrnd(intervalGamma);
			y 	   = binornd(n_true, observAlpha(iAlpha,:));

			maxObserv = max(y);
			iThreshOfMaxObserv = find(n_max >= ceil(1/observAlpha(iAlpha,1) * maxObserv), 1);

			%test pgffa
			tStart = tic;
			ll_pgffa_temp = gf_forward(y, intervalGamma, observAlpha(iAlpha,1), intervalDelta);
			llPGFFA(iAlpha)      = llPGFFA(iAlpha) + ll_pgffa_temp;
			runtimePGFFA(iAlpha) = runtimePGFFA(iAlpha) + toc(tStart);

			%test fa with each setting of n_max
			llEqualityMatched = false;
			for iThresh = 1:nThresh
				tStartPGFFA = tic;
				[~,~,ll_fa_temp] = forward_messages(y, intervalGamma, intervalDelta, observAlpha(iAlpha, :), n_max(iThresh));
                % ll_fa_temp = exp(ll_fa_temp);
				llFA(iAlpha, iThresh)      = llFA(iAlpha, iThresh) + ll_fa_temp;
				runtimeFAiter = toc(tStartPGFFA);
				runtimeFA(iAlpha, iThresh) = runtimeFA(iAlpha, iThresh) + runtimeFAiter;
				if iThresh == iThreshOfMaxObserv
					runtimeFAatMaxYByAlpha(iAlpha) = runtimeFAatMaxYByAlpha(iAlpha) + runtimeFAiter;
                end
				if ~llEqualityMatched && abs(ll_fa_temp - ll_pgffa_temp) <= llEpsilon
					runtimeFAatLLEquality(iAlpha) = runtimeFAatLLEquality(iAlpha) + runtimeFAiter;
					llEqualityMatched = true;
				end
			end
			if ~llEqualityMatched
				%if this happens rarely, it is probably a rounding error, or a numeric precision issue
				%try increasing epsilon if it is a problem
				warning('Warning: likelihood did not converge, using runtime of last iteration')
				runtimeFAatLLEquality(iAlpha) = runtimeFAatLLEquality(iAlpha) + runtimeFAiter;
			end

			runtimeIter = toc(tStart);
			runtimeTot  = sum(runtimeFA(:)) + sum(runtimePGFFA(:));
			runtimeMean = runtimeTot / ((iAlpha - 1) * nIter + iter);
			runtimeRem  = runtimeMean * ((nAlpha - iAlpha) * nIter + nIter - iter);

			% keyboard
			if ~silent
				if runtimeIter < 60
					fprintf(': %.1fs (~%.2fm remaining)\n', runtimeIter, runtimeRem/60);
				else
					fprintf(': %.2fm (~%.2fm remaining)\n', runtimeIter/60, runtimeRem/60);
				end
			end
		end
	end

	%average out all the runtimes and likelihoods by the number of iterations they were repeated
	runtimeFA              = runtimeFA              ./ nIter;
	runtimePGFFA           = runtimePGFFA           ./ nIter;
	runtimeFAatMaxYByAlpha = runtimeFAatMaxYByAlpha ./ nIter;
	runtimeFAatLLEquality  = runtimeFAatLLEquality  ./ nIter;
	llFA                   = llFA                   ./ nIter;
	llPGFFA                = llPGFFA                ./ nIter;

	iThreshOfTwiceNHat = find(n_max >= ceil(1.5*N_hat), 1);
	if isempty(iThreshOfTwiceNHat)
		iThreshOfTwiceNHat = nThresh;
	end
	runtimeFAatTwiceNHat   = runtimeFA(:, iThreshOfTwiceNHat);

	LINEWIDTH = 4;
	LABELFONT = 20;
	TITLEFONT = 22;
	LEGENDFONT = 12;
	if ~silent
		figure
		hold on

		plot(observAlpha(:,1), runtimeFAatTwiceNHat, 'LineWidth', LINEWIDTH)
		plot(observAlpha(:,1), runtimeFAatMaxYByAlpha, 'LineWidth', LINEWIDTH)
		plot(observAlpha(:,1), runtimeFAatLLEquality, 'LineWidth', LINEWIDTH)
		plot(observAlpha(:,1), runtimePGFFA, 'LineWidth', LINEWIDTH)

		title('\alpha vs runtime of FA and PGFFA', 'FontSize', TITLEFONT)
		xlabel('\alpha',                           'FontSize', LABELFONT)
		ylabel('mean runtime (s)',                 'FontSize', LABELFONT)

		legend({'$$n_{max} = 2\hat{N}$$'; ...
			    '$$n_{max} = \alpha^{-1}max(y)$$'; ...
			    '$$LL(FA) = LL(PGFFA)$$'; ...
			    '$$PGFFA$$'}, ...
			    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','East')
		hold off

		figure
		hold on

		plot(observAlpha(:,1), runtimeFAatMaxYByAlpha, 'LineWidth', LINEWIDTH)
		plot(observAlpha(:,1), runtimeFAatLLEquality, 'LineWidth', LINEWIDTH)
		plot(observAlpha(:,1), runtimePGFFA, 'LineWidth', LINEWIDTH)

		title('\alpha vs runtime of FA and PGFFA', 'FontSize', TITLEFONT)
		xlabel('\alpha',                           'FontSize', LABELFONT)
		ylabel('mean runtime (s)',                 'FontSize', LABELFONT)

		legend({'$$n_{max} = \alpha^{-1}max(y)$$'; ...
			    '$$LL(FA) = LL(PGFFA)$$'; ...
			    '$$PGFFA$$'}, ...
			    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','East')
		hold off

		figure
		hold on

		plot(observAlpha(:,1), runtimePGFFA, 'LineWidth', LINEWIDTH)

		title('\alpha vs runtime of PGFFA', 'FontSize', TITLEFONT)
		xlabel('\alpha',                    'FontSize', LABELFONT)
		ylabel('mean runtime (s)',          'FontSize', LABELFONT)

		% legend({'$$n_{max} = \alpha^{-1}max(y)$$'; ...
		% 	    '$$LL(FA) = LL(PGFFA)$$'; ...
		% 	    '$$PGFFA$$'}, ...
		% 	    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','East')
		hold off
	end
end