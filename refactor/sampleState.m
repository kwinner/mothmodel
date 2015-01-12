function [ y, state ] = sampleState( theta, params )
%SAMPLESTATE Generate a sample of transect counts given dist'n params
%   Inputs:
%       N      - number of individuals in population
%       mu     - emergence distribution mean
%       sigma  - emergence distribution stddev
%       lambda - lifespan distribution parameter
%       alpha  - observability parameter
%       t      - vector of observation times

%T = # of observations
T = length(params.t);

%compute p
state.p = ppdf(theta, params);

%% sample the q (outcome) values ~ Multinomial(N,p)
% state.q = sample_hist(state.p(:), params.N);
state.q = mnrnd(params.N, state.p(:));
state.q = reshape(state.q, T+1, T+1);

%% compute n (the abundance at each observation t]
state.n = abundancy(state.q);

%% sample each y_k (the number of individuals observed at each t)

y = binornd(state.n,params.alpha);

end

