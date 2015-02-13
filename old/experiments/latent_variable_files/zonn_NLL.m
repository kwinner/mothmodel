function NLL = zonn_NLL( y, T, mu, sigma, lambda, N, alpha, varargin )
%compute the NLL according to the classic zonneveld model

parser = inputParser;
addOptional  (parser, 'p_t',      []);

parser.parse(varargin{:});
p_t      = parser.Results.p_t;

%compute p_t, the prob that an indiv is alive at the times in T
if isempty(p_t)
	p_t = presence_probs(mu, sigma, lambda, T, inf);
end

NLL = sum(observation_NLL(y, N, p_t, alpha));

end

