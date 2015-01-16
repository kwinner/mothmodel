function NLL = latentvar_NLL( Q, T, mu, sigma, lambda, N, alpha, varargin )
%compute the NLL of the latent variables

parser = inputParser;
addOptional(parser, 'p_ij', []);

parse(parser, varargin{:});
p_ij = parser.Results.p_ij;

%compute p_ij if it wasn't already provided
if isempty(p_ij)
	p_ij = birthdeath_pmf(mu, sigma, lambda, N, T);
end

NLL = -sum(logmnpdf(Q(:), p_ij(:)));

end

