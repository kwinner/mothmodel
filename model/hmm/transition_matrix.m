function PT = transition_matrix( gamma_k, delta_k, n_max )
% TRANSITION_MATRIX := compute probabilities of transitioning from t_{k-1} to t_k
% PT = transition_matrix( gamma_k, delta_k, n_max )
%
% INPUTS
% required:
%    gamma_k = the mean arrival for new individuals from t_{k-1} to t_k
%    delta_k = the probability of an individual surviving from t_{k-1} to t_k
%    n_max   = maximum possible abundance at time t_k (positive int)
%
% OUTPUTS
%    PT      = transition matrix for all n_{k-1} to n_k

	assert(n_max > 0, 'n_max must be a positive int, was %s', n_max)

	PT = arrayfun(@(n_prev) transition_distn(gamma_k, delta_k, n_max, n_prev), ...
		          0:n_max, ...
		          'UniformOutput', false);
	PT = vertcat(PT{:});
end

