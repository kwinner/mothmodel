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
%		F'(1) = e^(a+b) * [f'(1) + af(1)]
f_p = f(2:end) .* (1:D);
F_p = exp(a+b) * (sum(f_p) + a * sum(f));
mean = F_p;

if nargout > 1
	%variance is given by:
    %       Var = F''(1) - [F'(1)]^2 + F'(1)
	%		F''(1) = e^(a+b) * [f''(1) + 2af'(1) + a^2 f(1)]
    f_pp = f_p(2:end) .* (1:D-1);
    F_pp = exp(a+b) * (sum(f_pp) + 2 * a * sum(f_p) + (a^2) * sum(f));
    var  = F_pp - (F_p ^ 2) + F_p;
end

end

