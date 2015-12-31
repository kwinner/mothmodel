function [psi, z, loglikelihood] = forward_messages( rateFunc, serviceDistn, T, y, N_hat, alpha, varargin )
% FORWARD_MESSAGES := compute the log-space normalized messages, psi.
% psi = presence_prob( arrivalDistn, serviceDistn, T, y, N_hat, alpha, varargin )
%      if only one output is used, the messages are restored before being returned
% [psi, z, likelihood] = presence_prob( arrivalDistn, serviceDistn, T, y, N_hat, alpha, varargin )
%      if > one output, the normalized messages are restored
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector [1 x K] of observation times (sample times)
%    y            = vector [1 x K] of observed counts
%    N_hat        = mean super population size (positive int)
%    alpha        = detection probability (probability)
% paramvalue:
%    'n_max'      = maximum possible abundance at each observation (positive int)
%
% OUTPUTS
%    psi           = matrix [n_max x K] of messages (either normalized or restored)
%                    messages will be normalized if z is returned, restored if not
%    z             = vector [1 x K] of normalizing constants for psi
%    loglikelihood = log of the joint likelihood of y

parser = inputParser;
addParamValue(parser, 'n_max',      ceil(N_hat))

parse(parser, varargin{:})
n_max      = parser.Results.n_max;

%for readability
K = numel(T);

psi = zeros(n_max+1, K);

%fill out the observation distribution
P_y = arrayfun(@(y_k) ...
	           binopdf(y_k, 0:n_max, alpha)' ...
	           , y ...
	           , 'UniformOutput', false);
P_y = horzcat(P_y{:});

%compute the initial message (with no transition matrix)
P_0 = initial_distn(rateFunc, serviceDistn, T(1), N_hat, 'n_max', n_max);
psi(:,1) = P_y(:,1) .* P_0';

%normalize the first messages
z(1)     = log(sum(psi(:,1)));
psi(:,1) = exp(-z(1)) .* psi(:,1);

% ll = eval_gf(exp(z(1)) .* psi(:,1)', (0:.2:1));
% figure
% plot(ll);


%compute the rest of the messages
for k = 2:K
	%compute the transition matrix
	PT = transition_matrix(rateFunc, serviceDistn, T(k-1), T(k), N_hat, 'n_max', n_max);
	
	%compute the message
	psi(:,k) = P_y(:,k) .* sum(PT .* repmat(psi(:,k-1), 1, n_max + 1), 1)';

	%normalize the message
	z(k)     = log(sum(psi(:,k)));
	psi(:,k) = exp(-z(k)) .* psi(:,k);

	% ll = eval_gf(exp(z(k)) .* psi(:,k)', (0:.2:1));
	% figure
	% plot(ll);
	
end

%if z isn't used in the output, de-normalize the output messages
if nargout == 1
	for i = 1:K
		psi(:,i) = exp(sum(z(1:i))) .* psi(:,i);
	end
%compute the likelihood if it was requested
elseif nargout == 3
	%the likelihood is the sum of all the final messages at t=K
	%this needs to be denormalized still
	%since that sum will always be 1 due to normalization,
	%the likelihood is just the denormalization constant for the final messages
	loglikelihood = sum(z);
end

end

