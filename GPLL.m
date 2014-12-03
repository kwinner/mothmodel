function LL = GPLL(theta, y, params)

%compute the matrix p (a multinomial with compact parameters theta)
p = ppdf(theta, params);

%for each observation, compute the probability of being alive at that time
%note, this is mathematically equivalent to finding abundance on p instead of q
T = numel(y);

%from abundancy.mat, but real abundancy needs to cut out rounding error
%each element of n is a 2D sum of a submatrix of q
pt = arrayfun(@(k) sum(sum(p(1:k, k+1:T+1))), 1:T);

%compute the expected number of individuals observed at each obs
obsMean = params.alpha * params.N .* pt;

%compute the covariance of the observations
obsCov = GPcov(pt, theta, params);

%compute the LL
LL = log(mvnpdf(y, obsMean, obsCov));

end