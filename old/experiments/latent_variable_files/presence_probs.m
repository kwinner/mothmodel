function [p_t] = presence_probs( mu, sigma, lambda, T, discretization )
%computes the probability that an individual will be alive at each time t in T

if ~exist('discretization','var')
	discretization = inf;
end

%K = number of samples
K = numel(T);

if discretization == inf
	%use numeric integration instead of convolution
	p_t = arrayfun(@(t) quadgk(@(s) normpdf(s, mu, sigma) .* (1 - expcdf(t-s, lambda)), -inf, t), T);
else
	%X is the points to compute f_s and S_z at
	X = min(T) - std(T) : K/discretization : max(T) + std(T);

	f_s = normpdf(X, mu, sigma);
	S_z = 1 - expcdf(X, lambda);

	p_t = conv(f_s, S_z);
end

end

