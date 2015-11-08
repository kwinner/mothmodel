function P_0 = initial_distn( rateFunc, serviceDistn, t_0, N_hat, varargin )
% INITIAL_DISTN := compute the initial distribution of abundance at t_0
% P_0 = initial_distn( rateFunc, serviceDistn, t_0, N_hat, varargin )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    t_0          = the first observation time (real scalar)
%    N_hat        = mean super population size (positive int)
% paramvalue:
%    'n_max'      = maximum possible abundance at each observation (positive int)
%    'normalized' = should the probabilities be normalized to 1?
%                   this can be used if n_max is too low for the probabilities of high n to approach 0
%                   but really, n_max should be increased if this is expected
%
% OUTPUTS
%    P_0          = probability distribution over possible abundance (from 0 to n_max)

parser = inputParser;
addParamValue(parser, 'n_max',      N_hat)
addParamValue(parser, 'normalized', false)

parse(parser, varargin{:})
n_max      = parser.Results.n_max;
normalized = parser.Results.normalized;

%mean of initial arrivals (who survive to t_0)
arrival_mean = N_hat.*quadgk(@(s) ...
	                         rateFunc(s) ...
	                         .* (1-serviceDistn.cdf(t_0 - s)) ...
	                         , -inf, t_0);
P_0 = poisspdf(0:n_max, arrival_mean);

%if n_max is high enough, normalization is not needed
if normalized
	P_0 = P_0 ./ sum(P_0);
end

end

