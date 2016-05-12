function [likelihood,a,b,f,messages] = gf_forward( y, gamma, alpha, delta )
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
% f(s)e^(as+b), where f(s) is polynomial
K = length(gamma);

messages(K) = struct('a', 0, 'b', 0, 'f', 0);

a = 0;
b = 0;
f = [1]; %a vector representation of the coefficients of f(u)
         %f is "little endian": the leftmost entry is the constant
         %i.e. f(k) is the coefficient of u^{k-1}

for k = 1:K
	[a, b]  = immigrants(a, b, gamma(k));
	[a, f]  = observed(y(k), a, f, alpha);
    
    [b, f]  = normalize(b, f);
    
    messages(k).a = a;
    messages(k).b = b;
    messages(k).f = f;
    
	%last iteration does not have survivors, it just adds immigrants and conditions on observ
	if k <= length(delta)
		[a, b, f] = survivors(a, b, f, delta(k));
	end
end

%evaluate the GF at 1
likelihood = sum(f) * exp(a + b);

end

function [bPrime, fPrime] = normalize(b, f) 

    const = max(f);
    fPrime = f / const;
    bPrime = b + log(const);

end


function [aPrime, bPrime] = immigrants(a, b, gamma_k) 
% IMMIGRANTS := incorporate the GF of immigrants between t_{k-1} and t_k

aPrime = a + gamma_k;
bPrime = b - gamma_k;

end

function [ result ] = compose_poly_horner_special( f, g )
% COMPOSE_POLY_HORNER_SPECIAL Compose a polynomial f with a linear function
% g using Horner's method.
n = numel(f);
result = f(n);
for i = n-1:-1:1
    result = [result*g(1) 0] + [f(i) result*g(2)];
end
end

function [ result ] = compose_poly_horner( f, g )
% COMPOSE_POLY_HORNER Compose two polynomials using Horner's method
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

function [aPrime, bPrime, fPrime] = survivors(a, b, f, delta_k)
% SURVIVORS := thin the population by delta_k

aPrime = a * delta_k;
bPrime = b + a * (1 - delta_k);

fPrime = compose_poly_horner_special(f, [1-delta_k, delta_k]);

end


function [aPrime, fPrime] = observed(y_k, a, f, alpha)
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
    denom = gamma(x+1)*gamma(y_k-x+1);    % x! (y_k - x) !
    
	% move (1-alpha) inside f() to the coefficients
    numCo = (1-alpha) .^ (0:(length(xthDerivative)-1));
    
    % multiply by a^(y_k - x)
    acPowers = a^(y_k - x);
    
	% add to running sum
	fPrime = left_aligned_sum(fPrime, xthDerivative .* numCo ./ denom * acPowers);
end

% multiply by alpha^y_k
fPrime = fPrime * alpha^y_k; %mult by (alpha)^y_k

% multiply by s^y_k
fPrime = [zeros(1,y_k), fPrime(:)']; %increase all powers by y_k

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

