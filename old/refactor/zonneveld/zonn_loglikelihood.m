function LL = zonn_loglikelihood( y, theta_emergence, theta_lifespan, t, alpha, N )
%ZONN_LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here

T = numel(t);
mu     = theta_emergence(1);
sigma  = theta_emergence(2);
lambda = theta_lifespan;

p_t = zeros(1,T);

for i = 1:T
	p_t(i) = quadgk(@(s) integrand(s, t(i), mu, sigma, lambda), -inf, t(i));
end

n_t = N .* p_t;

LL = sum(log(poisspdf(y, n_t.* alpha)));

end

function value = integrand(s, t, mu, sigma, lambda)

value = normpdf(s, mu, sigma) .* (1 - expcdf(t-s, lambda));

end