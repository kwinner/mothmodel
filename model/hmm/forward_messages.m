function [psi, z, loglikelihood] = forward_messages( y, gamma, delta, alpha, n_max )
% FORWARD_MESSAGES := compute the log-space normalized messages, psi
% psi = forward_messages( y, gamma, delta, alpha, n_max )
%     if only one output is used, the messages are restored before being returned
% [psi, z, loglikelihood] = forward_messages( y, gamma, delta, alpha, n_max )
%     if > 1 output, the normalized messages are restored
%
% INPUTS
% required:
%    y     = vector [1 x K]   of observed counts
%    gamma = vector [1 x K+1] of mean arrival rate between each observation
%            note: if gamma is a scalar, it is expanded to match y
%    delta = vector [1 x K] of survival probabilities between each observation
%            note: if delta is a scalar, it is expanded to match y
%    alpha = vector [1 x K] of observation probabilities at each observation
%            note: if alpha is a scalar, it is expanded to match y
%    n_max = maximum possible abundance at each observation (positive int)
%
% OUTPUTS
%    psi          = matrix [n_max x K] of messages (either normalized or restored)
%                   messages will be normalized if z is returned, restored if not
%    z            = vector [1 x K] of normalizing constants for psi
%    logliklihood = log of the joint likelihood of y given gamma, delta, alpha

	assert(n_max > 0, 'n_max must be a positive int, was %s', n_max)

	%expand gamma, delta, alpha if necessary
	K = numel(y); %K = number of observations
	if numel(gamma) == 1; gamma = repmat(gamma, 1, K + 1); end
	if numel(delta) == 1; delta = repmat(delta, 1, K); end
	if numel(alpha) == 1; alpha = repmat(alpha, 1, K); end

	%initialize messages
	psi = zeros(n_max + 1, K);

	%fill out the observation distribution
	P_y = arrayfun(@(y_k, alpha_k) ...
		           binopdf(y_k, 0:n_max, alpha_k)' ...
		           , y, alpha ...
		           , 'UniformOutput', false);
	P_y = horzcat(P_y{:});

	%compute the initial message (transition distn w/ no n_prev and no survivors)
	P_0      = transition_distn(gamma(1), 1, n_max, 0);
	psi(:,1) = P_y(:,1) .* P_0';

	%normalize the first messages
	z(1)     = log(sum(psi(:,1)));
	psi(:,1) = exp(-z(1)) .* psi(:,1);

	%compute the rest of the messages
	for k = 2:K
		%compute the transition matrix
		PT = transition_matrix(gamma(k), delta(k-1), n_max);

		%complete the message
		psi(:,k) = P_y(:,k) .* sum(PT .* repmat(psi(:,k-1), 1, n_max + 1), 1)';

		%normalize the message
		z(k)     = log(sum(psi(:,k)));
		psi(:,k) = exp(-z(k)) .* psi(:,k);
	end

	%if z isn't used in the output, de-normalize the output messages
	if nargout == 1
		for i = 1:K
			psi(:,i) = exp(sum(z(1:i))) .* psi(:,i);
		end
	%compute the likelihood if it was requested
	elseif nargout == 3
		%the likelihood is the sum of all the final messages at t=K
		%but since those messages still need to be denormalized,
		%the loglikelihood is just the sum of the denormalization constants
		loglikelihood = sum(z);
	end
end

