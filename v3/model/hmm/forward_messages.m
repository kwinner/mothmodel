function psi = forward_messages( rateFunc, serviceDistn, T, y, N_hat, alpha, varargin )

parser = inputParser;
addParamValue(parser, 'n_max',      N_hat)

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

%compute the rest of the messages
for k = 2:K
	%compute the transition matrix
	PT = transition_matrix(rateFunc, serviceDistn, T(k-1), T(K), N_hat, 'n_max', n_max);
	psi(:,k) = P_y(:,k) .* (PT * psi(:,k-1));
end

end

