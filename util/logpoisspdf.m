function y = logpoisspdf( x, lambda )
% LOGPOISSPDF := compute the poisson pdf at x w/ mean lambda, in log space
% y = logpoisspdf( x, lambda )
%
% INPUTS
% required:
%    x      = point to evaluate at
%    lambda = mean of the poisson
%
% OUTPUTS
%    y      = probability of x given lambda

y = log(lambda) .* x - gammaln(x+1) - lambda;

end

