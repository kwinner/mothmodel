function [tOptim, nIters, tMarg] = salamander()
	TICKLABELFONT = 14;
	LABELFONT  = 20;
	TITLEFONT  = 22;
	LEGENDFONT = 16;
	FIGSIZE    = [0 0 10 6];
	% MEANCOLOR  = pink(100); MEANCOLOR = MEANCOLOR(90, :);
	MEANCOLOR = [0.35, 0.35, 0.35];

	salamanderData = csvread('~/Work/Data/salamander.csv', 1, 1);

	nonzerosites = sum(salamanderData, 2) > 0;

	salamanderData = salamanderData(nonzerosites, :);

	T = [5 6 17 18 29 30 41 42 53 54 65 66 77 78];
	R = size(salamanderData, 1);
	K = size(salamanderData, 2);

	t = tic;

	optimAlg      = 'interior-point';
	options       = optimoptions('fmincon', 'Algorithm', optimAlg, 'Display', 'off');
	problem.objective = @(theta) objective1(theta, salamanderData, T);
	problem.x0        = [0.5, max(salamanderData(:)), 0.5];
	problem.lb        = [0, 1, 0];
	problem.ub        = [1, inf, inf];
	problem.solver    = 'fmincon';
	problem.options   = options;

	[theta_hat, ~, ~, min_output] = fmincon(problem);

	nIters = min_output.funcCount;
 
	alpha_hat      = theta_hat(1);
	N_hat          = theta_hat(2);
	survivProb_hat = theta_hat(3);

	tOptim = toc(t);

	rateFunc = makedist('uniform', 'lower', min(T)-12, 'upper', max(T)+12);
	rateFunc = @rateFunc.pdf;
	serviceDistn = makedist('exponential', 'mu', survivProb_hat); 
	gamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
	delta = survival_prob(serviceDistn, T);

	tMarg = zeros(R, K, 3);

	for r = 1:R
	% for r = 1
		y = salamanderData(r, :);

		[~, ~, ~, ~, messages] = gf_forward(y, gamma, alpha_hat, delta);

		for i = 1:K
			t = tic;
			[lll, a, b, f] = gf_tail_eliminate(y, gamma, alpha_hat, delta, i, messages(i).a, messages(i).b, messages(i).f);
			f = f';
			tMarg(r,i,1) = toc(t);

			t = tic;
			pmf(:,i) = pgf2pmf(f, a, b-lll, 'K', 40)';
			tMarg(r,i,2) = toc(t);

			t = tic;
			[meanVec(i), varVec(i)] = moments_pgf(f, a, b - lll);
			tMarg(r,i,3) = toc(t);
		end

		% pmf = flipud(pmf);

		h = figure('Visible', 'off', 'Units', 'inches', 'Position', FIGSIZE);
		% h = figure('Units', 'inches', 'Position', FIGSIZE);
		hold on
		colormap(pink(255));
		imagesc(pmf);
		scatter(1:K, salamanderData(r,:)./alpha_hat, 200, 'd', 'filled', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2)
		scatter(1:K, meanVec, 150, '+', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', MEANCOLOR, 'LineWidth', 4)
		% for i  = 1:K
		% 	plot([0.55+(i-1), 1.45+(i-1)], [meanVec(i), meanVec(i)], 'Color', [1 1 1], 'LineWidth', 5)
		% end
		% scatter(1:K, meanVec, 200, 'd', 'filled', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2)
		% [l, icons, plots] = legend({'Observed Counts'}, 'Position', [0.6428    0.7657    0.1988    0.0370], 'FontSize', LEGENDFONT);
		[l, icons, plots] = legend({'Observed Counts / $$\hat{\rho}$$', 'Mean'}, 'Position', [0.5797    0.8045    0.2649    0.0993], 'FontSize', LEGENDFONT, 'Interpreter', 'Latex');
		title(sprintf('Posterior marginals for site #%d', max(find(nonzerosites, r))), 'FontSize', TITLEFONT)
		xlabel('T', 'FontSize', LABELFONT)
		set(gca, 'TickDir', 'out')
		xlim([0.5, K + 0.5]);
		ylim([0.5, 40.5]);
		set(gca, 'YDir', 'Normal')
		set(gca, 'FontSize', TICKLABELFONT)
		set(gca, 'YTick', 0:5:40)
		set(gca, 'YTickLabel', 0:5:40)
		set(gca, 'XTick', 1:K)
		set(gca, 'XTickLabel', T)

		% keyboard

		print(h, sprintf('marginals_site%d', max(find(nonzerosites, r))), '-dpng');
		close(h)
	end

end

function nll = objective1(theta, y, T)
	R = size(y, 1);

	alpha = theta(1);
	N = theta(2);
	survivProb = theta(3);

	rateFunc = makedist('uniform', 'lower', min(T)-12, 'upper', max(T)+12);
	rateFunc = @rateFunc.pdf;
	serviceDistn = makedist('exponential', 'mu', survivProb);

	gamma = immigration_rate(rateFunc, serviceDistn, T, N);
	delta = survival_prob(serviceDistn, T);

	nll = 0;

	for r = 1:R
		nll = nll - gf_forward(y(r,:), gamma, alpha, delta);
	end

end

function nll = objective2(theta, y, T)
	R = size(y, 1);

	alpha = theta(1);
	N = theta(2);
	survivProb = theta(3);

	gamma = N .* ones(size(T));
	delta = survivProb .* ones(size(T));

	nll = 0;

	for r = 1:R
		nll = nll - gf_forward(y(r,:), gamma, alpha, delta);
	end

end

