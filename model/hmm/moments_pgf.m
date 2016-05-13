function [mean, var] = moments_pgf(f, a, b)
% MOMENTS_PGF := Compute the first two central moments of a PGF
% [mean, var] = moments_pgf(f, a, b)
%
% Compute the first two central moments (mean and variance) 
%  of a PGF of the form:
% 		F(s) = f(s)*exp(as+b)
% where f(s) is a polynomial in s of degree D
%
% INPUTS
% required:
%    f = a [1 x D+1] vector of polynomial coefficients
%        such that f[1] is the coefficient of s^0,
%        f[2] is the coefficient of s^1, etc.
%    a = a scalar
%    b = a scalar

D = numel(f) - 1;

%mean is just the first derivative of F(s) evaluated at s = 1
%which is given by:
%		F'(1) = e^(a+b) * sum[k=0:D]((a+k)f[k])
mean = exp(a+b) * sum(((a:D+a) .* f));

if nargout > 1
	%variance is the second derivative of F(s), evaluated at s = 1
	%which is given by:
	%		F''(s) = e^(a+b) * sum[k=0:D]{((a+k)^2-k)f[k]}
	var = exp(a+b) * sum((((a:D+a).^2-(0:D)) .* f))
end

end

