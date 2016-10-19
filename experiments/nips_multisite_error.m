%written 5/17/16
function [errorPGFFA, errorTrunc, coveragePGFFA, coverageTrunc, CI_widthPGFFA, CI_widthTrunc, alpha_hatPGFFA, N_hatPGFFA, alpha_hatTrunc, N_hatTrunc, divergePGFFA, divergeTrunc] = nips_multisite_error(varargin)
	NalphaProduct = 10;
	n_max         = 50;
	K             = 10;
    R             = 20;

	nIter         = 1;
	% alphaVec      = 0.1:0.2:1.00;
    alphaVec      = 0.05:0.05:.25;
	nAlpha        = numel(alphaVec);

	optimAlg      = 'interior-point';
	options       = optimoptions('fmincon', 'Algorithm', optimAlg, 'Display', 'off', 'StepTolerance', 1e-15);

	%prepare output
	alpha_hatPGFFA = zeros(nAlpha, nIter);
	N_hatPGFFA     = zeros(nAlpha, nIter);
	errorPGFFA     = zeros(nAlpha, 2, nIter);
	CI_widthPGFFA  = zeros(nAlpha, 2, nIter);
	coveragePGFFA  = zeros(nAlpha, 2, nIter);
	divergePGFFA   = zeros(nAlpha, nIter);
	alpha_hatTrunc = zeros(nAlpha, nIter);
	N_hatTrunc     = zeros(nAlpha, nIter);
	errorTrunc     = zeros(nAlpha, 2, nIter);
	CI_widthTrunc  = zeros(nAlpha, 2, nIter);
	coverageTrunc  = zeros(nAlpha, 2, nIter);
	divergeTrunc   = zeros(nAlpha, nIter);

	%now iterate over values of alpha
	for iAlpha = 1:nAlpha
		%split N, alpha to maintain NalphaProduct
		alpha = alphaVec(iAlpha);
		N     = round(NalphaProduct / alpha);

		%repeat nIter times
		for iter = 1:nIter
			%sample K observations
            n = poissrnd(N, R, 1);
			y = binornd(repmat(n, 1, K), alpha);

			while any(isnan(y(:))) || any(y(:) > n_max)
				warning('Generated observations which are impossible under n_max.')

            	n = poissrnd(N, R, 1);
				y = binornd(repmat(n, 1, K), alpha);				
			end

			fprintf('N = %d, alpha = %0.2f, y = %s\n', N, alpha, mat2str(y));

			%learn parameter estimates w/ pgffa
			pgffa.objective = @(theta) pgffa_objective(theta, y);
			pgffa.x0        = [0.5, max(y(:))*2];
			pgffa.lb        = [0, 1];
			pgffa.ub        = [1, inf];
			pgffa.solver    = 'fmincon';
			pgffa.options   = options;

			[theta_hat, ~, exitflag, ~, ~, ~, theta_hessian] = fmincon(pgffa);
			alpha_hatPGFFA(iAlpha, iter) = theta_hat(1);
			N_hatPGFFA(iAlpha, iter)     = theta_hat(2);

			if rcond(theta_hessian) < 1e-15
				divergePGFFA(iAlpha, iter) = 1;
				fprintf('PGFFA, exit w/ %d, infinite estimate\n', exitflag);
			else
				fprintf('PGFFA, exit w/ %d\n', exitflag);
			end

			errorPGFFA(iAlpha, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthPGFFA(iAlpha, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coveragePGFFA(iAlpha, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthPGFFA(iAlpha, :, iter);

			%learn parameter estimates w/ the truncated n mixture model
			trunc.objective = @(theta) trunc_objective(theta, y, n_max);
			trunc.x0        = [0.5, max(y(:))*2];
			trunc.lb        = [0, 1];
			trunc.ub        = [1, inf];
			trunc.solver    = 'fmincon';
			trunc.options   = options;

			[theta_hat, ~, exitflag, ~, ~, ~, theta_hessian] = fmincon(trunc);
			alpha_hatTrunc(iAlpha, iter) = theta_hat(1);
			N_hatTrunc(iAlpha, iter)     = theta_hat(2);

			if rcond(theta_hessian) < 1e-15
				divergeTrunc(iAlpha, iter) = 1;
				fprintf('Trunc, exit w/ %d, infinite estimate\n', exitflag);
			else
				fprintf('Trunc, exit w/ %d\n', exitflag);
			end

			errorTrunc(iAlpha, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthTrunc(iAlpha, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coverageTrunc(iAlpha, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthTrunc(iAlpha, :, iter);
		end
	end
end

function nll = pgffa_objective(theta, y)
    R = size(y, 1);
	K = size(y, 2);

	alpha = theta(1);
	N     = theta(2);

	gamma = [N, zeros(1, K-1)];
	delta = ones(K-1, 1);

    nll = zeros(1,R);
    for r = 1:R
        [nll(r), ~, ~, ~, ~] = gf_forward(y(r,:), gamma, alpha, delta);
    end
    nll = -sum(nll);
end

function nll = trunc_objective(theta, y, n_max)
    R = size(y, 1);

	nll = zeros(1,R);
    for r = 1:R
        for n = 0:n_max
            nll(r) = nll(r) + poisspdf(n, theta(2)) * prod(binopdf(y(r,:), n, theta(1)));
        end
    end
    nll = -sum(log(nll));

    if isinf(nll) || isnan(nll)
    	keyboard
    end
end

