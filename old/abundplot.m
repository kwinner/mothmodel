function abundplot(theta, params, y, t)

if exist('t', 'var')
	t_y = params.t;
	params.t = t;
end

p = ppdf(theta, params);

T = numel(params.t);

pt = arrayfun(@(k) sum(sum(p(1:k, k+1:T+1))), 1:T);
pt = params.alpha .* pt;

plot(params.t, params.N .* pt);

% pdf = @(s,t) normpdf(s, theta.mu, theta.sigma) .* exp(-(t-s).*theta.lambda);

% for t = 1:numel(params.t)
% 	pt(t) = params.N * params.alpha .* quadgk(@(s) pdf(s,params.t(t)), -inf, params.t(t));
% end

% plot(params.t, pt)

if exist('y', 'var')
	hold on
	plot(t_y, y, 'o')
end

end