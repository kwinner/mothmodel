function PT_n = transition_distn( rateFunc, serviceDistn, t_1, t_2, N_hat, n_1, varargin )
% TRANSITION_DISTN := compute probabilities of transitioning from n_1 individuals to n = 0->n_max indivs at time t_2
% PT_n = transition_distn( rateFunc, serviceDistn, t_1, t_2, N_hat, n_1, varargin )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    t_1          = time of conditioned observation (real scalar)
%    t_2          = time of predicted observation (real scalar)
%    N_hat        = mean super population size (positive int)
%    n_1          = number of individuals at time t_1 (nonnegative int)
% paramvalue:
%    'n_max'      = maximum possible abundance at each observation (positive int)
%    'normalized' = should the probabilities be normalized to 1?
%                   this can be used if n_max is too low for the probabilities of high n to approach 0
%                   but really, n_max should be increased if this is expected
%
% OUTPUTS
%    PT_n         = probability distribution over next states (for n_2 in 0:n_max)

parser = inputParser;
addParamValue(parser, 'n_max',      N_hat)
addParamValue(parser, 'normalized', false)

parse(parser, varargin{:})
n_max      = parser.Results.n_max;
normalized = parser.Results.normalized;

%compute the probability of [0, ..., N_max] individuals arriving in this interval
arrival_mean = N_hat .* quadgk(@(s) ...
	                           rateFunc(s) ...
	                           .* (1-serviceDistn.cdf(t_2 - s)) ...
	                           , t_1, t_2);
arrival_distn = poisspdf(0:n_max, arrival_mean);

%compute the probability of [0, ..., n_1] survivors from the previous observation
survival_distn = binopdf(0:n_1, n_1, (1-serviceDistn.cdf(t_2 - t_1)));

PT_n = conv(arrival_distn, survival_distn);

%trim the end of P_k
PT_n = PT_n(1:n_max+1);

%if n_max is high enough, normalization is not needed
if normalized
	PT_n = PT_n ./ sum(PT_n);
end

end

