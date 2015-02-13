function LL = zonn_LL( y, T, mu, sigma, lambda, N, alpha, varargin )
%compute the LL according to the classic zonneveld model

parser = inputParser;
addOptional  (parser, 'p_t',      []);

parser.parse(varargin{:});
p_t      = parser.Results.p_t;

%compute p_t, the prob that an indiv is alive at the times in T
if isempty(p_t)
	p_t = presence_probs(mu, sigma, lambda, T, inf);
end

LL = sum(observation_LL(y, N, p_t, alpha));

end

