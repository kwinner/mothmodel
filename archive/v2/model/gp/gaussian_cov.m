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

