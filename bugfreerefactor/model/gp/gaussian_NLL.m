function NLL = gaussian_NLL( y, T, mu, sigma, lambda, N, alpha, varargin )
%compute the NLL according to the classic zonneveld model

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
	NLL = -inf;
else
	NLL = -mvnormpdfln(y', obsMean', [], obsCov);
end

end


function Sigma = gaussian_cov( p_t, lambda, T, N, alpha )

K = numel(T);

Sigma = zeros(K,K);

%compute the variance (the diagonal of Sigma)
Sigma(logical(eye(K))) = N * (alpha .* p_t - alpha^2 .* p_t.^2);

%compute the covariance (the off diagonal entries)
for i = 1:K-1
	for j = i+1:K
		Sigma(i,j) = N * alpha^2 * p_t(i) * (exp(-(T(j) - T(i))/lambda) - p_t(j));
		Sigma(j,i) = Sigma(i,j);
	end
end

end

