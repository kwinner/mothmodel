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
	[a, b]       = immigrants(a, b, gamma(k));
	[a, f]       = observed(y(k), a, c, f, alpha);

	%last iteration does not have survivors, it just adds immigrants and conditions on observ
	if k <= length(delta)
		% [a, b, c, d] = survivors(a, b, c, d, delta(k));
		[a, b, c, d, f] = survivors2(a, b, c, d, f, delta(k));
	end

	% ll = arrayfun(@(s) sum(f .* ((c * s + d) .^ (0:length(f)-1))) .* exp(a * s + b), (0:.2:1));
	% figure
	% plot(ll);
	
end

%evaluate the GF at 1
likelihood = sum(f .* ((c + d) .^ (0:length(f)-1))) .* exp(a + b);

end


function [aPrime, bPrime] = immigrants(a, b, gamma_k) 
% IMMIGRANTS := incorporate the GF of immigrants between t_{k-1} and t_k

aPrime = a + gamma_k;
bPrime = b - gamma_k;

end


function [aPrime, bPrime, cPrime, dPrime] = survivors(a, b, c, d, delta_k)
% SURVIVORS := thin the population by delta_k

aPrime = a * delta_k;
bPrime = b + a * (1 - delta_k);
cPrime = c * delta_k;
dPrime = d + c * (1 - delta_k);

end

function [aPrime, bPrime, cPrime, dPrime, fPrime] = survivors2(a, b, c, d, f, delta_k)
% SURVIVORS := thin the population by delta_k

aPrime = a * delta_k;
bPrime = b + a * (1 - delta_k);

cPrime = c * delta_k;
dPrime = d + c * (1 - delta_k);

fPrime = compose_poly(f, [dPrime cPrime]);
cPrime = 1;
dPrime = 0;

end

function h = compose_poly(f, g)
% tic
f = poly2sym(fliplr(f));
g = poly2sym(fliplr(g));
h = compose(f, g);
h = fliplr(sym2poly(h));
% toc
end


function [aPrime, fPrime] = observed(y_k, a, c, f, alpha)
% OBSERVED := observe y_k, a value thinned by alpha
%note, these comments won't make much sense w/o the general form available in the paper

aPrime = a * (1 - alpha);

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
	denom = (a^x) * factorial(x) * factorial(y_k-x);

	%move (1-alpha) inside f() to the coefficients
	numCo = c^x .* (1-alpha)^x .* (1-alpha) .^ (0:(length(xthDerivative)-1));

	%add to running sum
	fPrime = left_aligned_sum(fPrime, xthDerivative .* numCo ./ denom);
end

%multiply the sum by (a * alpha)^y_k * s^y_k
fPrime = fPrime .* ((a * alpha) ^ y_k); %mult by (a*alpha)^y_k
fPrime = [zeros(1,y_k), fPrime]; %increase all powers by y_k

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

