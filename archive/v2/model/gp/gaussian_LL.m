function LL = gaussian_LL( y, T, mu, sigma, lambda, N, alpha, varargin )
%compute the LL according to the classic zonneveld model

parser = inputParser;
addOptional(parser, 'p_t', []);

parse(parser, varargin{:});
p_t      = parser.Results.p_t;

%compute p_t, the prob that an indiv is alive at the times in T
if isempty(p_t)
	p_t = presence_probs(mu, sigma, lambda, T, inf);
end

%convert p_t to the mean # of indiv alive at each t
obsMean = alpha * N .* p_t;

%compute the covariance of the observations
K = numel(T);
obsCov = gaussian_cov(p_t, lambda, T, N, alpha);

[~,p] = chol(obsCov);
if p ~= 0
	LL = -inf;
else
	LL = mvnormpdfln(y', obsMean', [], obsCov);
end

end

