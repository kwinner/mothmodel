function [y, S, Z, Q] = samplestate( mu, sigma, lambda, N, alpha, T )

%K = number of observations
K = numel(T);

%sample N birthtimes from Normal(mu, sigma)
S = normrnd(mu, sigma, N, 1);

%sample N lifespans from Exp(lambda)
Z = exprnd(lambda, N, 1);

%compute death times
D = S + Z;

y = arrayfun(@(t) sum(S <= t & D >= t), T);

if nargout >= 4
	%append -inf, inf to T for the interval computations
	T = [-inf, T, inf];
	[birth_intervals, death_intervals] = meshgrid(1:K+1,1:K+1);
	Q = arrayfun(@(i,j) sum(S >= T(i) & S <= T(i+1) & D >= T(j) & D <= T(j+1)), birth_intervals', death_intervals');
end

end

