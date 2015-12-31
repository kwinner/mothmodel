function PT_n = transition_distn( gamma_k, delta_k, n_max, n_prev )
% TRANSITION_DISTN := compute probability of transitioning from n_prev to all n in [0,n_max]
% PT_n = transition_distn( gamma_k, delta_k, n_max, n_prev )
%
% INPUTS
% required:
%    gamma_k = the mean arrival for new individuals from t_{k-1} to t_k
%    delta_k = the probability of an individual surviving from t_{k-1} to t_k
%    n_max   = maximum possible abundance at time t_k (positive int)
%    n_prev  = the abundance at time t_{k-1}
%
% OUTPUTS
%    PT_n    = probability distribution over all n_k (unnormalized)
% note: this function is meant to be called over all possible n_prev to fill out the full matrix

	%build the arrival and survival distributions from gamma_k, delta_k
	%note: arrivals are assumed poisson and survival is assumed binomial
	arrival_distn  = poisspdf(0:n_max,  gamma_k);
	survival_distn = binopdf (0:n_prev, n_prev, delta_k);

	%convolve the arrival and survival probabilities
	PT_n = conv(arrival_distn, survival_distn);

	%trim off PT_n above n_max
	%note: PT_n(1) corresponds to p(n_k = 0), so PT_n is length n_max + 1
	PT_n = PT_n(1:n_max+1);
end

