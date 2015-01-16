function [p_ij] = birthdeath_pmf( mu, sigma, lambda, N, T )
%computes the probability that an individual will be born in interval i and die in j

%K = number of samples
K = numel(T);

%append infs to the end of the sampling times
T = [-inf, T, inf];

%initialize p
p_ij = zeros(K+1);

for i = 1:K+1
	for j = i:K+1
		p_ij(i,j) = quadgk(@(s) normpdf(s, mu, sigma) .* (expcdf(T(j+1) - s, lambda) - expcdf(T(j) - s, lambda)), T(i), T(i+1));
		%keep values in the upper triangular from going to 0
		%NOTE: this is really just a shortcut to making the computation above more numerically stable
		if p_ij(i,j) == 0
			p_ij(i,j) = realmin;
		end
	end
end

end

