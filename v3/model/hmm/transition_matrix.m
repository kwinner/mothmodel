function PT = transition_matrix( rateFunc, serviceDistn, t_1, t_2, N_hat, varargin )

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

