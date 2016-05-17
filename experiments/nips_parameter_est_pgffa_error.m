%written 5/17/16
function [errorPGFFA, errorTrunc, coveragePGFFA, coverageTrunc, CI_widthPGFFA, CI_widthTrunc, alpha_hatPGFFA, N_hatPGFFA, alpha_hatTrunc, N_hatTrunc] = nips_parameter_est_pgffa_error(varargin)
	NalphaProduct = 10;
	n_max         = 100;
	K             = 3;

	nIter         = 5;
	alphaVec      = 0.01:0.01:1.00;
	nAlpha        = numel(alphaVec

	N_LIMIT       = 50000; %limit for the optimizer
	optimAlg      = 'interior-point';
	options       = optimoptions('fmincon', 'Algorithm', optimAlg, 'Display', 'iter');

	%prepare output
	alpha_hatPGFFA = zeros(nAlpha, nIter);
	N_hatPGFFA     = zeros(nAlpha, nIter);
	errorPGFFA     = zeros(nAlpha, 2, nIter);
	CI_widthPGFFA  = zeros(nAlpha, 2, nIter);
	coveragePGFFA  = zeros(nAlpha, 2, nIter);
	alpha_hatTrunc = zeros(nAlpha, nIter);
	N_hatTrunc     = zeros(nAlpha, nIter);
	errorTrunc     = zeros(nAlpha, 2, nIter);
	CI_widthTrunc  = zeros(nAlpha, 2, nIter);
	coverageTrunc  = zeros(nAlpha, 2, nIter);

	%now iterate over values of alpha
	for iAlpha = alphaVec
		%split N, alpha to maintain NalphaProduct
		alpha = alphaVec(iAlpha);
		N     = NalphaProduct / alpha;

		%repeat nIter times
		for iter = 1:nIter
			%sample K observations
			y = binornd(N, alpha, K, 1);

			%learn parameter estimates w/ pgffa
			pgffa.objective = @(theta) pgffa_objective(theta, y);
			pgffa.x0        = [alpha, N];
			pgffa.lb        = [0, max(y)];
			pgffa.ub        = [1, N_LIMIT];
			pgffa.solver    = 'fmincon';
			pgffa.options   = options;

			[theta_hat, ~, ~, ~, ~, ~, theta_hessian] = fmincon(pgffa);
			alpha_hatPGFFA(iAlpha, iter) = theta_hat(1);
			N_hatPGFFA(iAlpha, iter)     = theta_hat(2);

			errorPGFFA(iAlpha, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthPGFFA(iAlpha, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coveragePGFFA(iAlpha, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthPGFFA(iAlpha, :, iter);

			%learn parameter estimates w/ the truncated n mixture model
			trunc.objective = @(theta) trunc_objective(theta, y);
			trunc.x0        = [alpha, N];
			trunc.lb        = [0, max(y)];
			trunc.ub        = [1, N_LIMIT];
			trunc.solver    = 'fmincon';
			trunc.options   = options;

			[theta_hat, ~, ~, ~, ~, ~, theta_hessian] = fmincon(trunc);
			alpha_hatTrunc(iAlpha, iter) = theta_hat(1);
			N_hatTrunc(iAlpha, iter)     = theta_hat(2);

			errorTrunc(iAlpha, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthTrunc(iAlpha, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coverageTrunc(iAlpha, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthTrunc(iAlpha, :, iter);
		end
	end
end

function nll = pgffa_objective(theta, y)
	K = numel(y);

	alpha = theta(1);
	N     = theta(2);

	gamma = [N, zeros(K-1, 1)];
	delta = ones(K-1, 1);

	nll = -gf_backward(y, gamma, alpha, delta);
end

function nll = trunc_objective(theta, y)
	nll = 0;
	for n = 0:n_max
		nll = nll - poisspdf(n, theta(2)) * prod(binopdf(y, n, theta(1)));
	end
end

