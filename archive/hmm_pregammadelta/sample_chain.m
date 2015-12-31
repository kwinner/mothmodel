function n = sample_chain( rateFunc, serviceDistn, T, N_hat, varargin )
% SAMPLE_CHAIN := draw one sample of abundance from the markov chain
% n = sample_chain( rateFunc, serviceDistn, T, N_hat )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector [1 x K] of observation times (sample times)
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
addParamValue(parser, 'n_max',      N_hat)

parse(parser, varargin{:})
n_max      = parser.Results.n_max;

n = zeros(size(T));

%sample the # of initial individuals
P_0 = initial_distn(rateFunc, serviceDistn, T(1), N_hat, 'n_max', n_max);
n(1) = randsample(0:n_max,1,true,P_0);

%compute transitions
for k = 2:numel(T)
	tk1 = T(k-1);
	tk2 = T(k);

	%sample transition
	PT_k = transition_distn(rateFunc, serviceDistn, tk1, tk2, N_hat, n(k-1), 'n_max', n_max);
	n(k) = randsample(0:n_max,1,true,PT_k);
end

end

