function [ y, q, n, p ] = sampleState( N, mu, sigma, lambda, alpha, t )
%SAMPLESTATE Generate a sample of transect counts given dist'n params
%   Inputs:
%       N      - number of individuals in population
%       mu     - emergence distribution mean
%       sigma  - emergence distribution stddev
%       lambda - lifespan distribution parameter
%       alpha  - observability parameter
%       t      - vector of observation times

%T = # of observations
T = length(t);

%initialize returns
y = zeros(1,T);
q = zeros(T+1,T+1);
n = zeros(1,T);
p = zeros(T+1,T+1);

%create emergence and lifespan dist'n objects
epdf = @(S) normpdf(S,mu,sigma);
lcdf = @(Z) expcdf(Z,lambda);
lifespanProb = @(a, b) max(0, lcdf(b) - lcdf(a));

%% compute p (by integrating the emer. and lifespan dist'ns)

%see doc for p distribution
%p is a matrix of outcome probabilities for the different birth and
% death intervals

intervalLimits = [-Inf, t, Inf];

for i=1:T+1
    for j=1:T+1
        %limits of emergence
        a1 = intervalLimits(i);
        b1 = intervalLimits(i+1);
        
        %limits of lifespan
        a2 = intervalLimits(j);
        b2 = intervalLimits(j+1);
        
        integrand = @(s) epdf(s).*lifespanProb(a2-s, b2-s);
        
        p(i,j) = quadgk(integrand, a1, b1);
    end
end

%renormalize p
p(:) = p(:) ./ sum(p(:));

%% sample the q (outcome) values ~ Multinomial(N,p)

q = mnrnd(N,p(:));
q = reshape(q, T+1, T+1);

%% compute n (the abundance at each observation t]

n = abundancy(q);

%% sample each y_k (the number of individuals observed at each t)

y = binornd(n,alpha);

end

