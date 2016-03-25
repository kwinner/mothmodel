function [beta, z, loglikelihood] = backward_messages( y, gamma, delta, alpha, n_max )
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
	beta = zeros(n_max + 1, K); %NOTE A ZEROTH INDEX IS ADDED TO BETA LATER
	beta(:,K) = ones(n_max + 1,1);

	%fill out the observation distribution
	P_y = arrayfun(@(y_k, alpha_k) ...
		           binopdf(y_k, 0:n_max, alpha_k)' ...
		           , y, alpha ...
		           , 'UniformOutput', false);
	P_y = horzcat(P_y{:});

	%normalize the initial messages
	z(K)     = log(sum(beta(:,K)));
	beta(:,K) = exp(-z(K)) .* beta(:,K);

	%compute the rest of the messages
	for k = K-1:-1:1
		%compute the transition matrix
		PT = transition_matrix(gamma(k+1), delta(k), n_max);

		%complete the message
		% beta(:,k) = P_y(:,k+1) .* sum(PT .* repmat(psi(:,k-1), 1, n_max + 1), 1)';
		beta(:,k) = PT * (beta(:,k+1) .* P_y(:,k+1));

		%normalize the message
		z(k)     = log(sum(beta(:,k)));
		beta(:,k) = exp(-z(k)) .* beta(:,k);
	end

	%add the zeroth message
	beta = [zeros(n_max + 1, 1), beta];
	PT = transition_distn(gamma(1), 1, n_max, 0);
	beta(:,1) = PT' .* beta(:,2) .* P_y(:,1);
	%normalize the last message
	z         = [log(sum(beta(:,1))), z];
	beta(:,1) = exp(-z(1)) .* beta(:,1);

	%if z isn't used in the output, de-normalize the output messages
	if nargout == 1
		for i = 1:K
			beta(:,i) = exp(sum(z(i:end))) .* beta(:,i);
		end
	%compute the likelihood if it was requested
	elseif nargout == 3
		%the likelihood is the sum of all the final messages at t=K
		%but since those messages still need to be denormalized,
		%the loglikelihood is just the sum of the denormalization constants
		loglikelihood = sum(z);
	end
end

