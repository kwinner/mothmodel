function y = logpoisspdf( x, lambda )

y = log(lambda) .* x - gammaln(x+1) - lambda;

end

