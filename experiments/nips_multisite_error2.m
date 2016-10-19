%written 5/17/16
function [errorPGFFA, errorTrunc, coveragePGFFA, coverageTrunc, CI_widthPGFFA, CI_widthTrunc, alpha_hatPGFFA, N_hatPGFFA, alpha_hatTrunc, N_hatTrunc, divergePGFFA, divergeTrunc, yRecord] = nips_multisite_error2(varargin)
	NalphaProduct = 10;
	n_max         = 50;
	K             = 10;
    R             = 20;

	nIter         = 100;
	lambdaVec     = 10:10:200;
	% lambdaVec     = 10:10:50;
	nLambda       = numel(lambdaVec);

	optimAlg      = 'interior-point';
	options       = optimoptions('fmincon', 'Algorithm', optimAlg, 'Display', 'off', 'StepTolerance', 1e-15);

	%prepare output
	alpha_hatPGFFA = zeros(nLambda, nIter);
	N_hatPGFFA     = zeros(nLambda, nIter);
	errorPGFFA     = zeros(nLambda, 2, nIter);
	CI_widthPGFFA  = zeros(nLambda, 2, nIter);
	coveragePGFFA  = zeros(nLambda, 2, nIter);
	divergePGFFA   = zeros(nLambda, nIter);
	alpha_hatTrunc = zeros(nLambda, nIter);
	N_hatTrunc     = zeros(nLambda, nIter);
	errorTrunc     = zeros(nLambda, 2, nIter);
	CI_widthTrunc  = zeros(nLambda, 2, nIter);
	coverageTrunc  = zeros(nLambda, 2, nIter);
	divergeTrunc   = zeros(nLambda, nIter);

	%now iterate over values of alpha
	for iLambda = 1:nLambda
		%split N, alpha to maintain NalphaProduct
		N = lambdaVec(iLambda);
		alpha = NalphaProduct / N;

		%repeat nIter times
		for iter = 1:nIter
			%sample K observations
            n = poissrnd(N, R, 1);
			y = binornd(repmat(n, 1, K), alpha);
			yRecord(iLambda, iter, :) = y(:);

			while any(isnan(y(:))) || any(y(:) > n_max)
				warning('Generated observations which are impossible under n_max.')

            	n = poissrnd(N, R, 1);
				y = binornd(repmat(n, 1, K), alpha);				
			end

			fprintf('N = %d, alpha = %0.2f, y = %s\n', N, alpha, mat2str(y));

			%learn parameter estimates w/ pgffa
			pgffa.objective = @(theta) pgffa_objective(theta, y);
			pgffa.x0        = [0.2, max(y(:))*2];
			pgffa.lb        = [0, 1];
			pgffa.ub        = [1, inf];
			pgffa.solver    = 'fmincon';
			pgffa.options   = options;

			[theta_hat, ~, exitflag, ~, ~, ~, theta_hessian] = fmincon(pgffa);
			alpha_hatPGFFA(iLambda, iter) = theta_hat(1);
			N_hatPGFFA(iLambda, iter)     = theta_hat(2);

			if rcond(theta_hessian) < 1e-15
				divergePGFFA(iLambda, iter) = 1;
				fprintf('PGFFA, exit w/ %d, infinite estimate\n', exitflag);
			else
				fprintf('PGFFA, exit w/ %d\n', exitflag);
			end

			errorPGFFA(iLambda, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthPGFFA(iLambda, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coveragePGFFA(iLambda, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthPGFFA(iLambda, :, iter);

			%learn parameter estimates w/ the truncated n mixture model
			trunc.objective = @(theta) trunc_objective(theta, y, n_max);
			trunc.x0        = [0.5, max(y(:))*2];
			trunc.lb        = [0, 1];
			trunc.ub        = [1, inf];
			trunc.solver    = 'fmincon';
			trunc.options   = options;

			[theta_hat, ~, exitflag, ~, ~, ~, theta_hessian] = fmincon(trunc);
			alpha_hatTrunc(iLambda, iter) = theta_hat(1);
			N_hatTrunc(iLambda, iter)     = theta_hat(2);

			if rcond(theta_hessian) < 1e-15
				divergeTrunc(iLambda, iter) = 1;
				fprintf('Trunc, exit w/ %d, infinite estimate\n', exitflag);
			else
				fprintf('Trunc, exit w/ %d\n', exitflag);
			end

			errorTrunc(iLambda, :, iter) = (theta_hat - [alpha, N]) .^ 2;

			theta_cov                      = inv(theta_hessian);
			CI_widthTrunc(iLambda, :, iter) = abs(2.*sqrt(diag(theta_cov)))';

			coverageTrunc(iLambda, :, iter) = abs(theta_hat - [alpha, N]) <= CI_widthTrunc(iLambda, :, iter);
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
    	temp = -inf(n_max - max(y(r,:)),1);
        for n = max(y(r,:)):n_max
            temp(n) = logpoisspdf(n, theta(2)) + sum(logbinopdf(y(r,:), n, theta(1)));
            % nll(r) = nll(r) + poisspdf(n, theta(2)) * prod(binopdf(y(r,:), n, theta(1)));
        end
        nll(r) = logsumexp(temp);
    end
    nll = -sum(nll);

    if isinf(nll) || isnan(nll)
    	keyboard
    end
end


function s = logsumexp(a, dim)
% Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
% Default is dim = 1 (columns).
% logsumexp(a, 2) will sum across rows instead of columns.
% Unlike matlab's "sum", it will not switch the summing direction
% if you provide a row vector.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

if nargin < 2
  dim = 1;
end

% subtract the largest in each column
[y, i] = max(a,[],dim);
dims = ones(1,ndims(a));
dims(dim) = size(a,dim);
a = a - repmat(y, dims);
s = y + log(sum(exp(a),dim));
i = find(~isfinite(y));
if ~isempty(i)
  s(i) = y(i);
end

end