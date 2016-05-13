%written 5/10/16
function [runtimePGFFA, runtimeFAatMaxYByAlpha, runtimeFAatLLEquality, runtimeFAatTwiceNHat] = fa_pgffa_runtime_vs_N(varargin)
	DEFAULT_NITER		= 25;
	DEFAULT_N_HAT       = [10:10:100, 125:25:500];
	% DEFAULT_N_HAT = [125:25:500];
	% DEFAULT_N_HAT       = 10:10:100;
	DEFAULT_N_LIMIT     = max(DEFAULT_N_HAT*1.5);

	%actually, n_max can't really be set this way, but maybe a parameter controlling step size?
	% DEFAULT_N_MAX		= 1:floor(1/DEFAULT_N_HAT*20):DEFAULT_N_HAT * 1.5 + 1;

	DEFAULT_K           = 0;
	DEFAULT_MU          = 8;
	DEFAULT_SIGMA       = 4;
	DEFAULT_LAMBDA      = 3;
	DEFAULT_T           = 1:4:20;
	DEFAULT_ALPHA       = 0.5;
	DEFAULT_GAMMA       = [];
	DEFAULT_DELTA       = [];
	DEFAULT_SILENT      = false;
	DEFAULT_LL_EPSILON  = 0.001;

	parser = inputParser;
	addParamValue(parser, 'nIter',     DEFAULT_NITER)
	addParamValue(parser, 'alpha',     DEFAULT_ALPHA)
	addParamValue(parser, 'K',         DEFAULT_K)
	addParamValue(parser, 'mu',        DEFAULT_MU)
	addParamValue(parser, 'sigma',     DEFAULT_SIGMA)
	addParamValue(parser, 'lambda',    DEFAULT_LAMBDA)
	addParamValue(parser, 'T',         DEFAULT_T)
	addParamValue(parser, 'N_hat',     DEFAULT_N_HAT)
	addParamValue(parser, 'n_max_lim', DEFAULT_N_LIMIT)
	% addParamValue(parser, 'n_max',   DEFAULT_N_MAX)
	addParamValue(parser, 'gamma',     DEFAULT_GAMMA)
	addParamValue(parser, 'delta',     DEFAULT_DELTA)
	addParamValue(parser, 'silent',    DEFAULT_SILENT)
	addParamValue(parser, 'epsilon',   DEFAULT_LL_EPSILON)

	parse(parser, varargin{:})
	nIter         = parser.Results.nIter;
	observAlpha   = parser.Results.alpha;
	K             = parser.Results.K;
	mu            = parser.Results.mu;
	sigma         = parser.Results.sigma;
	lambda        = parser.Results.lambda;
	T             = parser.Results.T;
	N_hat         = parser.Results.N_hat;
	n_max_lim	  = parser.Results.n_max_lim;
	% n_max         = parser.Results.n_max;
	unscaledGamma = parser.Results.gamma;
	intervalDelta = parser.Results.delta;
	silent        = parser.Results.silent;
	llEpsilon     = parser.Results.epsilon;

	%populate per-interval params
	if K == 0
		K = numel(T);
	end
	%fa supports per-sample values for alpha, but we are assuming a constant detection rate
	%so alpha needs to be expanded
	if size(observAlpha, 2) == 1
		observAlpha = repmat(observAlpha, 1, K);
	end
	if isempty(unscaledGamma == []) && isempty(intervalDelta == [])
		arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
		rateFunc     = @arrivalDistn.pdf;

		%service distribution
		serviceDistn = makedist('Exp', 'mu', lambda);

		unscaledGamma = immigration_rate(rateFunc, serviceDistn, T, 1);
		intervalDelta = survival_prob(serviceDistn, T);
	elseif isempty(unscaledGamma == []) || isempty(intervalDelta == [])
		error('gamma and delta must both be specified or unspecified.')
	else
		warning('Note: this experiment takes /unscaled/ gamma, delta parameters. Do not multiply by N_hat')
	end

	nN = numel(N_hat);

	%the complete FA results cannot be stored trivially, as the number of values of n_max that will be tested is variable
	% runtimeFA              = zeros(nN, nThresh);
	runtimePGFFA           = zeros(nN, 1);
	runtimeFAatMaxYByAlpha = zeros(nN, 1); %has to be tracked separately because corresponding iThresh depends on ymax
	runtimeFAatLLEquality  = zeros(nN, 1); %also would depend on y
	runtimeFAatTwiceNHat   = zeros(nN, 1);
	% llFA		           = zeros(nN, nThresh);
	llPGFFA                = zeros(nN, 1);
	runtimeTot			   = 0;
	for iN = 1:nN
		if ~silent, fprintf('---Experiment %d of %d, N_hat = %d---\n', iN, nN, N_hat(iN)); end

		%scale gamma by N(iN)
		intervalGamma = unscaledGamma .* N_hat(iN);

		for iter = 1:nIter
			if ~silent, fprintf('Iteration %d of %d', iter, nIter); end

			%sample observations
			n_true = poissrnd(intervalGamma);
			y 	   = binornd(n_true, observAlpha(1));

			n_max_of_max_observ = ceil(1/observAlpha(1) * max(y));
			n_max_of_2N         = ceil(N_hat(iN) * 1.5);

			%test pgffa
			tStart = tic;
			ll_pgffa_temp = gf_forward(y, intervalGamma, observAlpha(1), intervalDelta);
			llPGFFA(iN)      = llPGFFA(iN) + ll_pgffa_temp;
			runtimePGFFA(iN) = runtimePGFFA(iN) + toc(tStart);

			%test fa with each setting of n_max
			n_max = max(1,min(y)); %n_max is basically a loop iterator, but the increment is nonstandard
			n_max_step = 1;
			llEqualityMatched = false;
			while(n_max <= n_max_lim)
				tStartPGFFA = tic;
				[~,~,ll_fa_temp] = forward_messages(y, intervalGamma, intervalDelta, observAlpha, n_max);
                % ll_fa_temp = exp(ll_fa_temp);
				runtimeFAiter = toc(tStartPGFFA);
				if n_max == n_max_of_max_observ
					runtimeFAatMaxYByAlpha(iN) = runtimeFAatMaxYByAlpha(iN) + runtimeFAiter;
                end
                if n_max == n_max_of_2N
                	runtimeFAatTwiceNHat(iN)   = runtimeFAatTwiceNHat(iN) + runtimeFAiter;
                end
				if ~llEqualityMatched && abs(ll_fa_temp - ll_pgffa_temp) <= llEpsilon
					runtimeFAatLLEquality(iN) = runtimeFAatLLEquality(iN) + runtimeFAiter;
					llEqualityMatched = true;
				end
				%if all 3 of these are true, then we should have sampled all 3 of our strategies above
				if llEqualityMatched && n_max >= n_max_of_max_observ && n_max >= n_max_of_2N
					break
				else
					%convoluted n_max increment
					if ~llEqualityMatched
						%if the ll of FA still hasn't converged, then increment n_max normally
						n_max = n_max + n_max_step;
					elseif n_max < min(n_max_of_max_observ, n_max_of_2N)
						%otherwise, we can just jump to n_max_of_max_observ and n_max_of_2N
						%first, we'll jump to the lesser of the two
						n_max = min(n_max_of_max_observ, n_max_of_2N);
					else
						%then, jump to the second
						n_max = max(n_max_of_max_observ, n_max_of_2N);
						%note, the loop should execute one more time here with this last n_max, then exit above
					end
				end
			end
			if ~llEqualityMatched
				%if this happens rarely, it is probably a rounding error, or a numeric precision issue
				%try increasing epsilon if it is a problem
				%also, in this version, the only way for this to happen would be if n_max > n_max_limit, which would take... ages.
				warning('Warning: likelihood did not converge, using runtime of last iteration')
				runtimeFAatLLEquality(iN) = runtimeFAatLLEquality(iN) + runtimeFAiter;
			end

			runtimeIter = toc(tStart);
			runtimeTot  = runtimeTot + runtimeIter;
			runtimeMean = runtimeTot / ((iN - 1) * nIter + iter);
			runtimeRem  = runtimeMean * ((nN - iN) * nIter + nIter - iter);

			% keyboard
			%this "remaining time" printout is a lot less useful in this case, since n_max grows with each iteration
			%maybe if we did the different values of n_max in reverse order? start with the biggest?
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
	runtimePGFFA           = runtimePGFFA           ./ nIter;
	runtimeFAatMaxYByAlpha = runtimeFAatMaxYByAlpha ./ nIter;
	runtimeFAatLLEquality  = runtimeFAatLLEquality  ./ nIter;
	runtimeFAatTwiceNHat   = runtimeFAatTwiceNHat   ./ nIter;
	llPGFFA                = llPGFFA                ./ nIter;

	LINEWIDTH = 4;
	LABELFONT = 20;
	TITLEFONT = 22;
	LEGENDFONT = 12;
	if ~silent
		figure
		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(N_hat, runtimeFAatTwiceNHat, 'LineWidth', LINEWIDTH, 'Color', colOrd(4,:))
		plot(N_hat, runtimeFAatMaxYByAlpha, 'LineWidth', LINEWIDTH, 'Color', colOrd(3,:))
		plot(N_hat, runtimeFAatLLEquality, 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(N_hat, runtimePGFFA, 'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		title('$$\hat{N}$$ vs runtime of FA and PGFFA', 'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		xlabel('$$\hat{N}$$',                           'FontSize', LABELFONT, 'Interpreter', 'Latex')
		ylabel('mean runtime (s)',                      'FontSize', LABELFONT)

		legend({'$$n_{max} = 2\hat{N}$$'; ...
			    '$$n_{max} = \alpha^{-1}max(y)$$'; ...
			    '$$LL(FA) = LL(PGFFA)$$'; ...
			    '$$PGFFA$$'}, ...
			    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		figure
		hold on

		plot(N_hat, runtimeFAatMaxYByAlpha, 'LineWidth', LINEWIDTH, 'Color', colOrd(3,:))
		plot(N_hat, runtimeFAatLLEquality, 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(N_hat, runtimePGFFA, 'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		title('$$\hat{N}$$ vs runtime of FA and PGFFA', 'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		xlabel('$$\hat{N}$$',                           'FontSize', LABELFONT, 'Interpreter', 'Latex')
		ylabel('mean runtime (s)',                      'FontSize', LABELFONT)

		legend({'$$n_{max} = \alpha^{-1}max(y)$$'; ...
			    '$$LL(FA) = LL(PGFFA)$$'; ...
			    '$$PGFFA$$'}, ...
			    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		figure
		hold on

		plot(N_hat, runtimePGFFA, 'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		title('$$\hat{N}$$ vs runtime of PGFFA', 'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		xlabel('$$\hat{N}$$',                    'FontSize', LABELFONT, 'Interpreter', 'Latex')
		ylabel('mean runtime (s)',               'FontSize', LABELFONT)

		% legend({'$$n_{max} = \alpha^{-1}max(y)$$'; ...
		% 	    '$$LL(FA) = LL(PGFFA)$$'; ...
		% 	    '$$PGFFA$$'}, ...
		% 	    'Interpreter','Latex','FontSize',LEGENDFONT,'Location','East')
		hold off
	end
end