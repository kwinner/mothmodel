function PT = transition_matrix( rateFunc, serviceDistn, t_1, t_2, N_hat, varargin )
% TRANSITION_MATRIX := compute probabilities of transitioning from n_1 individuals to n = 0->n_max indivs at time t_2
% PT = transition_matrix( rateFunc, serviceDistn, t_1, t_2, N_hat, varargin )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    t_1          = time of conditioned observation (real scalar)
%    t_2          = time of predicted observation (real scalar)
%    N_hat        = mean super population size (positive int)
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
%normalized = parser.Results.normalized;

%compute all the rows of the transition matrix
PT = arrayfun(@(n_0) ...
	          transition_distn(rateFunc, serviceDistn, t_1, t_2, N_hat, n_0, varargin{:}) ...
	          , 0:n_max ...
	          , 'UniformOutput', false);

%concatenate the rows into the final transition matrix
PT = vertcat(PT{:});

end

