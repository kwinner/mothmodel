function P_0 = initial_distn( rateFunc, serviceDistn, t_0, N_hat, varargin )

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

