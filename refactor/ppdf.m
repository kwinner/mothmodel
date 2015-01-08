function [P] = ppdf(theta, params)

%factor theta
if isstruct(theta)
	mu = theta.mu;
	sigma = theta.sigma;
	lambda = theta.lambda;
else
	mu = theta(1);
	sigma = theta(2);
	lambda = theta(3);
end

%initialize P
T = numel(params.t);
P = zeros(T+1);

%append -inf and inf to t, to allow individuals to be born
%anytime before the first observ or after the last observ
params.t = [-inf, params.t, inf];

%loop over i,j: the upper triangular half of P
for i = 1:T+1
	for j = i:T+1
		%integrate s from t(i) to t(i+1)
		%s is time of birth
		P(i,j) = quadgk(@(s) integrand(i,j, s, mu, sigma, lambda, params), params.t(i), params.t(i+1));
	end
end

end

%the meat and potatoes, prob of birth * prob of lifespan
function [value] = integrand(i, j, s, mu, sigma, lambda, params)
DEFAULT_LOGISTIC_MIX = 0;

[zmin, zmax] = lifespan_domain(s, params.t(j), params.t(j+1));
value = ((1-DEFAULT_LOGISTIC_MIX).*normpdf(s, mu, sigma)+DEFAULT_LOGISTIC_MIX.*pdf('logistic',s,mu,sigma)) .* (expcdf(zmax, lambda) - expcdf(zmin, lambda));
end