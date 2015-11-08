function likelihood = dail_madsen(y, lambda, alpha, N_max)

likelihood = 0;
for n = min(y):N_max
	% likelihood = likelihood + poisspdf(n, lambda) .* prod(arrayfun(@(yi) binopdf(yi,n,alpha),y));
	likelihood = likelihood + poisspdf(n, lambda) .* prod(binopdf(y,n,alpha));
end

end

% function pdf = binopdf(x, n, p)
% 	pdf = exp(gammaln(n+1) - gammaln(x+1) - gammaln(n-x+1) + x.*log(p) + (n-x).*log(1-p));
% end