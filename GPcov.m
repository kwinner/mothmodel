function obsCov = GPcov (pt, theta, params)

T = numel(params.t);

%during optimization, theta is sometimes a struct and sometimes a vector, this handles it
if isnumeric(theta)
	lambda = theta(3);
else
	lambda = theta.lambda;
end

%compute the covariance of the observations
obsCov = zeros(T,T);
for i = 1:T
	for j = i:T
		if i ~= j
			obsCov(i,j) = params.N * params.alpha^2 * pt(i) * (exp(-1/lambda * (params.t(j) - params.t(i))) - pt(j));
			obsCov(j,i) = obsCov(i,j);
		else
			obsCov(i,i) = params.N * (params.alpha * pt(i) - params.alpha^2 * pt(i)^2);
		end
	end
end

end