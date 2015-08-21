function PT_n = transition_distn( rateFunc, serviceDistn, t_1, t_2, N_hat, n_0, varargin )

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

%compute the probability of [0, ..., n_0] survivors from the previous observation
survival_distn = binopdf(0:n_0, n_0, (1-serviceDistn.cdf(t_2 - t_1)));

PT_n = conv(arrival_distn, survival_distn);

%trim the end of P_k
PT_n = PT_n(1:n_max+1);

%if n_max is high enough, normalization is not needed
if normalized
	PT_n = PT_n ./ sum(PT_n);
end

end

