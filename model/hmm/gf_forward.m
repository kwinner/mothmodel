function [likelihood,a,b,c,d,f] = gf_forward( y, gamma, alpha, delta )
% GF_FORWARD := forward algorithm for counts implemented with generating functions
% likelihood = gf_forward( y, gamma, alpha, delta )
%
% INPUTS
% required:
%    y     = vector [1 x K] of observed counts
%    gamma = vector [1 x K] of arrival rate
%    alpha = detection probability
%    delta = vector [1 x K] of death probability
%
% OUTPUTS
%    likelihood = exact joint likelihood of y given theta = {gamma, alpha, delta}

%at each step, we assume the generating function maintains this form
% f(cs+d)e^(as+b), where f(u) is polynomial
a = 0;
b = 0;
c = 1;
d = 0;
f = [1]; %a vector representation of the coefficients of f(u)
         %f is "little endian": the leftmost entry is the constant
         %i.e. f(k) is the coefficient of u^{k-1}

for k = 1:length(gamma)
	[a, b]        = immigrants(a, b, gamma(k));
	[a, c, d, f]  = observed(y(k), a, c, d, f, alpha);
    
    [b, f] = normalize(b, f);
    
	%last iteration does not have survivors, it just adds immigrants and conditions on observ
	if k <= length(delta)
		[a, b, c, d, f] = survivors(a, b, c, d, f, delta(k));
	end
end

%evaluate the GF at 1
likelihood = sum(f .* ((c + d) .^ (0:length(f)-1))) * exp(a + b);

end

function [bPrime, fPrime] = normalize(b, f) 

    c = max(f);
    fPrime = f / c;
    bPrime = b + log(c);

end


function [aPrime, bPrime] = immigrants(a, b, gamma_k) 
% IMMIGRANTS := incorporate the GF of immigrants between t_{k-1} and t_k

aPrime = a + gamma_k;
bPrime = b - gamma_k;

end

function [ result ] = compose_poly_horner( f, g )
% COMPOSE_POLY_HORNER Compose two polynomials using Horner's method
%    
%   h = compose_poly_horner(f, g)
%
n = numel(f);
result = f(n);
for i = n-1:-1:1
    result = conv(result, g);
    result(1) = result(1) + f(i);
end
end

% Compute composition of two polynomials by converting them to symbolic
% functions and then back 
function h = compose_poly(f, g)
f = poly2sym(fliplr(f));
g = poly2sym(fliplr(g));
h = compose(f, g);
h = fliplr(sym2poly(h));
end

function [aPrime, bPrime, cPrime, dPrime, fPrime] = survivors(a, b, c, d, f, delta_k)
% SURVIVORS := thin the population by delta_k

aPrime = a * delta_k;
bPrime = b + a * (1 - delta_k);

cPrime = c * delta_k;
dPrime = d + c * (1 - delta_k);

fPrime = compose_poly_horner(f, [dPrime cPrime]);
cPrime = 1;
dPrime = 0;

end


function [aPrime, cPrime, dPrime, fPrime] = observed(y_k, a, c, d, f, alpha)
% OBSERVED := observe y_k, a value thinned by alpha
%note, these comments won't make much sense w/o the general form available in the paper

aPrime = a * (1 - alpha);
cPrime = c;
dPrime = d;

%sum over x from 0 to y_k the xth derivative of f over a^x * x! * (y_k - x)!
fPrime = 0;
xthDerivative = f;
for x = 0:min(y_k,length(f)-1)
	if x > 0
		%compute the next partial of xthDerivative
		xthDerivative = xthDerivative(2:end); %throw away the constant
		xthDerivative = xthDerivative .* (1:length(xthDerivative)); %mult by power
	end

	%compute denominator
	denom = factorial(x) * factorial(y_k-x);

	%move (1-alpha) inside f() to the coefficients
    %numCo = ones(size(xthDerivative));
    %  BUG: this only works if d = 0
    numCo = (1-alpha) .^ (0:(length(xthDerivative)-1));
    
    % multiply by a^(y_k - x) times c^x
    acPowers = a^(y_k - x)*c^x;
    
	%add to running sum
	fPrime = left_aligned_sum(fPrime, xthDerivative .* numCo ./ denom * acPowers);
end

% % Convert from polynomial in t = cs+d to polynomial in s
% s = linspace(0,1,length(fPrime))';
% A = fliplr(vander(s));
% t = c*s + d;
% B = fliplr(vander(t));
% v = B * fPrime';
% fPrime = (A\v)';
% cPrime = 1;
% dPrime = 0;

%multiply by alpha^y_k
fPrime = fPrime * alpha^y_k; %mult by (alpha)^y_k

%multiply by s^y_k
fPrime = [zeros(1,y_k), fPrime(:)']; %increase all powers by y_k

%cPrime = 1-alpha;
%  BUG: this only works if d = 0
%
%  Otherwise it looks like we need to either express s^k as a polynomial in
%  t for t=(cs+d), or multiply out the polynomial in t using Binomial
%  expansion.
% 

end

function [C] = left_aligned_sum(A,B)
% LEFTALIGNEDSUM := add two vectors of different lengths by padding the shorter vector with zeros

aLength = length(A);
bLength = length(B);

if aLength < bLength
	C = [A, zeros(1,bLength - aLength)] + B;
elseif bLength < aLength
	C = [B, zeros(1,aLength - bLength)] + A;
else
	C = A + B;
end

end


function gamma = immigration_rate( rateFunc, serviceDistn, T, N_hat )
% IMMIGRATION_RATE := evaluate the rate of immigration between each observation time
% gamma = immigration_rate( rateFunc, serviceDistn, T )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector of observation times (reals)
%    N_hat        = 
%
% OUTPUTS
%    gamma        = vector of rates of new "successful" immigrants

%rate of arrivals who survive to t_0
gamma = N_hat .* arrayfun(@(t1,t2) ...
	                               quadgk(@(s) ...
	                                      rateFunc(s) ...
	                                      .* (1-serviceDistn.cdf(t2 - s)) ...
	                                      , t1, t2) ...
	                      , [-inf, T(1:end-1)], T);

end


function delta = survival_prob( serviceDistn, T )
% SURVIVAL_PROB := evaluate the probability of an individual surviving from T(k) to T(k+1)
% delta = survival_prob( serviceDistn, T )
%
% INPUTS
% required:
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector of observation times (reals)
%
% OUTPUTS
%    delta        = vector of survival probability of existing indivs

%probability of an individual living from t1 to t2
delta = arrayfun(@(t1,t2) (1-serviceDistn.cdf(t2 - t1)) ...
	             , T(1:end-1), T(2:end));

end

