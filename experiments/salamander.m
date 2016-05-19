function theta_hat = salamander()
	salamanderData = csvread('~/Work/Data/salamander.csv', 1, 1);

	salamanderData = salamanderData(sum(salamanderData, 2) > 0, :);

	T = [5 6 17 18 29 30 41 42 53 54 65 66 77 78];
	R = size(salamanderData, 1);
	K = size(salamanderData, 2);

	optimAlg      = 'interior-point';
	options       = optimoptions('fmincon', 'Algorithm', optimAlg, 'Display', 'iter');
	problem.objective = @(theta) objective1(theta, salamanderData, T);
	problem.x0        = [0.5, max(salamanderData(:)), 0.5];
	problem.lb        = [0, 1, 0];
	problem.ub        = [1, inf, inf];
	problem.solver    = 'fmincon';
	problem.options   = options;

	theta_hat = fmincon(problem);
 
	alpha_hat      = theta_hat(1);
	N_hat          = theta_hat(2);
	survivProb_hat = theta_hat(3);

	rateFunc = makedist('uniform', 'lower', min(T)-12, 'upper', max(T)+12);
	rateFunc = @rateFunc.pdf;
	serviceDistn = makedist('exponential', 'mu', survivProb_hat); 
	gamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
	delta = survival_prob(serviceDistn, T);

	% for r = 1:R
	for r = 1
		y = salamanderData(r, :);

		[~, ~, ~, ~, messages] = gf_forward(y, gamma, alpha_hat, delta);

		for i = 1:K
			[~, a, b, f] = gf_tail_eliminate(y, gamma, alpha_hat, delta, i, messages(i).a, messages(i).b, messages(i).f);
			f = f';

			pmf(:,i) = pgf2pmf(f, a, b, 'K', 60)';
		end
	end

	imagesc(pmf);
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

