function [pmf, rem] = pgf2pmf(f, a, b, varargin)
% PGF2PMF := Extract a pmf from the corresponding pgf
% pmf        = pgf2pmf(f, a, b, ...)
% [pmf, rem] = pgf2pmf(f, a, b, ...)
%
% Compute the coefficients of the power series of a PGF F(s) with form:
% 		F(s) = f(s)*exp(as+b)
% where f(s) is a polynomial in s of degree D
% since it is impossible to compute the infinite number of terms of the pmf,
% a limit is required. There are two ways to specify this limit:
% 		epsilon: keep computing terms until the total probability mass
%                accounted for is >= 1 - epsilon
%		K:       compute at least the first K terms
% if both are specified, we will compute terms until both conditions are satisfied
% if neither are given, we will use a default value of epsilon, w/ no limit on K
%
% INPUTS
% required:
%    f          = a [1 x D+1] vector of polynomial coefficients
%                 such that f[1] is the coefficient of s^0,
%                 f[2] is the coefficient of s^1, etc.
%    a          = a scalar
%    b          = a scalar
% paramvalue:
%    epsilon    = [0,1) a limiting parameter on how much probability mass needs to be accounted for
%    K          = Z^+   a limiting parameter on how many terms need to be computed
%    normalized = bool  a flag for whether to return the normalized or unnormalized pmf
%    stop_on_0  = bool  a flag for whether or not to stop when the terms of the pmf go to 0
%
% OUTPUTS
%    pmf        = the pmf (normalized or unnormalized)
%                 note: if K is specified, size(pmf) will be at least [1 x K]
%                       if K spec and not eps, then size(pmf) will be exactly [1 x K]
%    rem        = the amount of unaccounted for probability mass
%                 theoretically equal to 1 - sum(pmf)

DEFAULT_EPSILON    = 0.0001;
DEFAULT_NORMALIZED = true;
DEFAULT_STOP_ON_0  = false;

parser = inputParser;
%note that epsilon and K initially default to 0
%if neither is specified, DEFAULT_EPSILON will be used
addParamValue(parser, 'epsilon',    0.)
addParamValue(parser, 'K',          0)
addParamValue(parser, 'normalized', DEFAULT_NORMALIZED)
addParamValue(parser, 'stop_on_0',  DEFAULT_STOP_ON_0);
parser.parse(varargin{:});
epsilon    = parser.Results.epsilon;
K          = parser.Results.K;
normalized = parser.Results.normalized;
stop_on_0  = parser.Results.stop_on_0;
if K == 0 && epsilon == 0, epsilon = DEFAULT_EPSILON; end

%D = degree of f
D = numel(f) - 1;

%init rem, pmf
rem = 1.;
pmf = zeros(1, K);

%each term in pmf depends on this function g(x) = f[x]./a^x
g = f ./ (a .^ (0:D));

k = 0;
while (K == 0 || k < K) && (epsilon == 0 || rem >= epsilon)
	%throughout the comments, we will treat vectors as 0-indexed, 
	%and indicate the distinction by using square brackets to index in the comments
	%pmf[k] = a^k * e^b * sum{x=0:k}(g(x)/(k-x)!)

	pmf(k+1) = sum(g ./ gamma(k - (0:D) + 1)); %sum(g(x)/(k-x)!)

	%mult by a^k*e^b
	pmf(k+1) = pmf(k+1) * a^k * exp(b);

	%increment and recompute rem    
	rem = rem - pmf(k+1);
	k = k + 1;

	if stop_on_0 && pmf(k) == 0
		break
    end
end

if normalized
	%normalize pmf before function returns
% 	if rem < 1 && rem > 0
% 		pmf = pmf ./ rem;
% 	else
		pmf = pmf ./ sum(pmf);
% 	end
end

end

