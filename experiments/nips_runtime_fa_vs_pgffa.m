%written 5/16/16
function [runtimePGFFA, runtimeFAatOracle, runtimeFAatPoiss, runtimeFAatNegBin, n_maxOracle, n_maxPoiss, n_maxNegBin] = nips_runtime_fa_vs_pgffa(varargin)
	DEFAULT_NITER		  = 25;     %how many repetitions of the experiments to do

	%experimental parameters (will be iterated over)
	%population size
	DEFAULT_N_HAT         = [10:10:100, 125:25:500];
	%detection probability
	DEFAULT_ALPHA         = [0.05:0.05:1];

	DEFAULT_MU            = 8;      %mean of arrivals
	DEFAULT_SIGMA         = 4;      %sd of arrivals
	DEFAULT_LAMBDA        = 3;      %mean lifespan
	DEFAULT_T             = 1:4:20; %sampling times
	DEFAULT_GAMMA         = [];     %arrival rates
	DEFAULT_DELTA         = [];     %survival rates
	                                %if gamma and delta are empty, they'll be set from the Zonneveld style mu/sigma/lambda
	DEFAULT_SILENT        = false;  %whether to print progress statements
	DEFAULT_LL_EPSILON    = 0.001;  %how far off likelihood has to be to consider it to have "converged"
	DEFAULT_TAIL_PROB_LIM = 1e-10;  %how much of the probability mass allowed to be in the tail of the dist'n

	username = getenv('USER');
	if username == 'kwinner'
		DEFAULT_RESULT_DIR_BASE = '~/Work/Data/Results';
	else
		DEFAULT_RESULT_DIR_BASE = 'Results';
	end
	DEFAULT_RESULT_DIR = fullfile(DEFAULT_RESULT_DIR_BASE, datestr(now, 'yyyymmddTHHMMSSFFF'));

	parser = inputParser;
	addParamValue(parser, 'nIter',         DEFAULT_NITER)
	addParamValue(parser, 'N_hat',         DEFAULT_N_HAT)
	addParamValue(parser, 'alpha',         DEFAULT_ALPHA)
	addParamValue(parser, 'mu',            DEFAULT_MU)
	addParamValue(parser, 'sigma',         DEFAULT_SIGMA)
	addParamValue(parser, 'lambda',        DEFAULT_LAMBDA)
	addParamValue(parser, 'T',             DEFAULT_T)
	addParamValue(parser, 'gamma',         DEFAULT_GAMMA)
	addParamValue(parser, 'delta',         DEFAULT_DELTA)
	addParamValue(parser, 'silent',        DEFAULT_SILENT)
	addParamValue(parser, 'll_epsilon',    DEFAULT_LL_EPSILON)
	addParamValue(parser, 'tail_prob_lim', DEFAULT_TAIL_PROB_LIM);
	addParamValue(parser, 'result_dir',    DEFAULT_RESULT_DIR);

	parse(parser, varargin{:})
	nIter         = parser.Results.nIter;
	N_hatVec      = parser.Results.N_hat;
	alphaVec      = parser.Results.alpha;
	mu            = parser.Results.mu;
	sigma         = parser.Results.sigma;
	lambda        = parser.Results.lambda;
	T             = parser.Results.T;
	unscaledGamma = parser.Results.gamma;
	delta         = parser.Results.delta;
	silent        = parser.Results.silent;
	llEpsilon     = parser.Results.ll_epsilon;
	tailProbLim   = parser.Results.tail_prob_lim;
	resultDir     = parser.Results.result_dir;

	%prepare output directory
	if exist(resultDir, 'dir') ~= 7
		try
			mkdir(resultDir)
		catch exception
			error('Error: result directory does not exist and cannot be created.')
		end
	end

	%populate per-interval params
	K = numel(T);

	if isempty(unscaledGamma == []) && isempty(delta == [])
		arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
		rateFunc     = @arrivalDistn.pdf;

		%service distribution
		serviceDistn = makedist('Exp', 'mu', lambda);

		unscaledGamma = immigration_rate(rateFunc, serviceDistn, T, 1);
		delta         = survival_prob(serviceDistn, T);
	elseif isempty(unscaledGamma == []) || isempty(delta == [])
		error('Error: gamma and delta must both be specified or unspecified.')
	else
		warning('this experiment takes /unscaled/ gamma, delta parameters. Do not multiply by N_hat')
	end

	nN     = numel(N_hatVec);
	nAlpha = numel(alphaVec);

	%mostly just used for output purposes
	nExperiments = nN * nAlpha;
	iExperiment  = 1;

	%result initialization
	runtimePGFFA      = zeros(nN, nAlpha, nIter);
	runtimeFAatOracle = zeros(nN, nAlpha, nIter);
	runtimeFAatPoiss  = zeros(nN, nAlpha, nIter);
	runtimeFAatNegBin = zeros(nN, nAlpha, nIter);
	n_maxOracle       = zeros(nN, nAlpha, nIter); %record the values of n_max used for oracle
	n_maxPoiss        = zeros(nN);                %record the values of n_max used for poisson (depends only on N_hat)
	n_maxNegBin       = zeros(nN, nAlpha, nIter); %record the values of n_max used for negbin
	runtimeTot        = 0;
	for iN = 1:nN
		N_hat = N_hatVec(iN);

		%scale gamma by N(iN)
		gamma = unscaledGamma .* N_hat;

		%compute the n_max for the poisson tail guarantee
		%this could probably be vectorized a lot better
		%first we compute the maximum rate given by the following pattern:
		% lambda_1 = gamma_1
		% lambda_2 = gamma_1*delta_1 + gamma_2
		% lambda_3 = gamma_1*delta_1*delta_2 + gamma_2*delta_2 + gamma_3
		n_maxPoiss(iN) = max(arrayfun(@(k) sum(arrayfun(@(i) gamma(i) * prod(delta(i:k-1)) ...
			                                            , 1:k)) ...
			                          , 1:K) ...
			            );
		%then we evaluate the inverse poisson cdf at this maximum rate
		n_maxPoiss(iN) = ceil(poissinv(1-tailProbLim, n_maxPoiss(iN)));

		for iAlpha = 1:nAlpha
			alpha = alphaVec(iAlpha);

			if ~silent, fprintf('---Experiment %d of %d, N_hat = %d, alpha = %0.2f---\n', iExperiment, nExperiments, N_hat, alpha); end

			for iter = 1:nIter
				if ~silent, fprintf('Iteration %d of %d', iter, nIter); end

				%sample observations
				n_true = poissrnd(gamma);
				y 	   = binornd(n_true, alpha);

				%compute the n_max for the nbin observation tail guarantee
				n_maxNegBin(iN, iAlpha, iter) = ceil(nbininv (1-tailProbLim, max(y), alpha));

				%test pgffa
				tStart = tic;
				llPGFFA                        = gf_forward(y, gamma, alpha, delta);
				runtimePGFFA(iN, iAlpha, iter) = toc(tStart);

				%test fa with each setting of n_max
				%n_max starts at the minimum possible value, then increments carefully until ll converges,
				%at which point it skips ahead to do n_maxPoiss and n_maxNegBin before completing
				n_max     = max(1,min(y));
				converged = false;
				while(n_max <= max(n_maxPoiss(iN), n_maxNegBin(iN, iAlpha, iter)))
					%note loop will generally terminate at the test at the bottom, but this will prevent it
					%from running forever if likelihood fails to converge (probably due to imprecision)
					tStartFA = tic;
					[~,~,llFA] = forward_messages(y, gamma, delta, alpha, n_max);
					runtimeFA = toc(tStartFA);

					%update the three records as appropriate
					if n_max == n_maxPoiss(iN)
						runtimeFAatPoiss(iN, iAlpha, iter)  = runtimeFA;
	                end
	                if n_max == n_maxNegBin(iN, iAlpha, iter)
	                	runtimeFAatNegBin(iN, iAlpha, iter) = runtimeFA;
	                end
					if ~converged && abs(llFA - llPGFFA) <= llEpsilon
						n_maxOracle(iN, iAlpha, iter)       = n_max;
						runtimeFAatOracle(iN, iAlpha, iter) = runtimeFA;
						converged = true;
					end

					%convoluted n_max increment
					if ~converged
						%if the ll of FA still hasn't converged, then increment n_max normally
						n_max = n_max + 1;
					elseif n_max < min(n_maxPoiss(iN), n_maxNegBin(iN, iAlpha, iter))
						%otherwise, we can just jump to n_max_of_max_observ and n_max_of_2N
						%first, we'll jump to the lesser of the two
						n_max = min(n_maxPoiss(iN), n_maxNegBin(iN, iAlpha, iter));
					elseif n_max < max(n_maxPoiss(iN), n_maxNegBin(iN, iAlpha, iter))
						%then, jump to the second
						n_max = max(n_maxPoiss(iN), n_maxNegBin(iN, iAlpha, iter));
						%note, the loop should execute one more time here with this last n_max, then exit above
					else
						%then we're done here.
						break
					end
				end
				if ~converged
					%if this happens rarely, it is probably a rounding error, or a numeric precision issue
					%try increasing epsilon if it is a problem
					%also, in this version, the only way for this to happen would be if n_max > n_max_limit, which would take... ages.
					warning('Warning: likelihood did not converge, using runtime of last iteration')
					runtimeFAatOracle(iN, iAlpha, iter) = runtimeFA;
				end

				runtimeIter = toc(tStart);
				runtimeTot  = runtimeTot  + runtimeIter;
				runtimeMean = runtimeTot  / ((iExperiment - 1) * nIter + iter);
				runtimeRem  = runtimeMean * ((nExperiments - iExperiment) * nIter + nIter - iter);

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
			
			iExperiment = iExperiment + 1;
		end
	end

	%take the average over all iterations (for plotting)
	meanRuntimePGFFA      = mean(runtimePGFFA,      3);
	meanRuntimeFAatOracle = mean(runtimeFAatOracle, 3);
	meanRuntimeFAatPoiss  = mean(runtimeFAatPoiss,  3);
	meanRuntimeFAatNegBin = mean(runtimeFAatNegBin, 3);

	%write a meta file with the experimental parameters
	metaFile = fopen(fullfile(resultDir, 'meta.txt'), 'w');
	fprintf(metaFile, 'nIter=%d\n', nIter);
	fprintf(metaFile, 'mu=%0.2f\n', mu);
	fprintf(metaFile, 'sigma=%0.2f\n', sigma);
	fprintf(metaFile, 'lambda=%0.2f\n', lambda);
	fprintf(metaFile, 'T=%s\n', mat2str(T));
	fprintf(metaFile, 'gamma=%s\n', mat2str(unscaledGamma));
	fprintf(metaFile, 'delta=%s\n', mat2str(delta));
	fprintf(metaFile, 'll_epsilon=%e\n', llEpsilon);
	fprintf(metaFile, 'tail_prob_lim=%e', tailProbLim);
	fclose(metaFile);

	LINEWIDTH  = 4;
	LABELFONT  = 20;
	TITLEFONT  = 22;
	LEGENDFONT = 12;
	FIGSIZE    = [0 0 8 4];
	%plot the alpha vs runtime plots (one for each val of N_hat)
	for iN = 1:nN
		N_hat = N_hatVec(iN);

		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(alphaVec, meanRuntimeFAatPoiss(iN, :),  'LineWidth', LINEWIDTH, 'Color', colOrd(4,:))
		plot(alphaVec, meanRuntimeFAatNegBin(iN, :), 'LineWidth', LINEWIDTH, 'Color', colOrd(3,:))
		plot(alphaVec, meanRuntimeFAatOracle(iN, :), 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\alpha$$ vs runtime of FA and PGFFA, $$\\hat{N} = %d$$', N_hat), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)
		title(sprintf('alpha vs runtime of FA and PGFFA, N = %d', N_hat), ...
		      'FontSize', TITLEFONT)
		xlabel('alpha',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'FA - Poiss'; ...
			    'FA - NegBin'; ...
			    'FA - Oracle'; ...
			    'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('AllMethods_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
		close(hFig)

		%repeat for oracle vs PGFFA only
		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(alphaVec, meanRuntimeFAatOracle(iN, :), 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\alpha$$ vs runtime of FA and PGFFA, $$\\hat{N} = %d$$', N_hat), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)
		title(sprintf('alpha vs runtime of FA and PGFFA, N = %d', N_hat), ...
		      'FontSize', TITLEFONT)
		xlabel('alpha',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'FA - Oracle'; ...
			    'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('Oracle_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
		close(hFig)

		%repeat for PGFFA only
		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\alpha$$ vs runtime of PGFFA, $$\\hat{N} = %d$$', N_hat), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)
		title(sprintf('alpha vs runtime of PGFFA, N = %d', N_hat), ...
		      'FontSize', TITLEFONT)
		xlabel('alpha',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('PGFFA_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
		close(hFig)
	end

	%plot the N_hat vs runtime plots (one for each val of alpha)
	for iAlpha = 1:nAlpha
		alpha = alphaVec(iAlpha);

		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(N_hatVec, meanRuntimeFAatPoiss(:, iAlpha),  'LineWidth', LINEWIDTH, 'Color', colOrd(4,:))
		plot(N_hatVec, meanRuntimeFAatNegBin(:, iAlpha), 'LineWidth', LINEWIDTH, 'Color', colOrd(3,:))
		plot(N_hatVec, meanRuntimeFAatOracle(:, iAlpha), 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\hat{N}$$ vs runtime of FA and PGFFA, $$\\alpha = %0.2f$$', alpha), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)
		title(sprintf('N vs runtime of FA and PGFFA, alpha = %d', alpha), ...
		      'FontSize', TITLEFONT)
		xlabel('N',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'FA - Poiss'; ...
			    'FA - NegBin'; ...
			    'FA - Oracle'; ...
			    'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('AllMethods_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
		close(hFig)

		%repeat for oracle vs pgffa
		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(N_hatVec, meanRuntimeFAatOracle(:, iAlpha), 'LineWidth', LINEWIDTH, 'Color', colOrd(2,:))
		plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\hat{N}$$ vs runtime of FA and PGFFA, $$\\alpha = %0.2f$$', alpha), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)
		title(sprintf('N vs runtime of FA and PGFFA, alpha = %d', alpha), ...
		      'FontSize', TITLEFONT)
		xlabel('N',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'FA - Oracle'; ...
			    'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('Oracle_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
		close(hFig)

		%repeat for pgffa only
		% hFig = figure('Visible', 'off');
		% hFig.PaperUnits = 'inches';
		% hFig.PaperPosition = FIGSIZE;
		hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);

		hold on
		colOrd = get(gca, 'ColorOrder');

		plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINEWIDTH, 'Color', colOrd(1,:))

		% title(sprintf('$$\\hat{N}$$ vs runtime of PGFFA, $$\\alpha = %0.2f$$', alpha), ...
		%       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
		% xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
		% ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		title(sprintf('N vs runtime of PGFFA, alpha = %d', alpha), ...
		      'FontSize', TITLEFONT)
		xlabel('N',       'FontSize', LABELFONT)
		ylabel('mean runtime (s)', 'FontSize', LABELFONT)

		legend({'PGFFA'}, ...
			    'FontSize',LEGENDFONT,'Location','Northwest')
		hold off

		print(hFig, ...
			  fullfile(resultDir, ...
			  	       sprintf('PGFFA_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
		close(hFig)
	end
end