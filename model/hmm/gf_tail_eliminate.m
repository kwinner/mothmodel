
function [loglikelihood, a,b,f] = gf_tail_eliminate( y, lambda, rho, delta, i, a, b, f)
% GF_TAIL_ELIMINATE := eliminate variables i+1 through K and return the
%                      unnormalized marginal over variable i, given an
%                      initial factor over variable i (the alpha messages)
%
% [a,b,f] = gf_tail_eliminate( y, lambda, rho, delta, i, a, b, f)
%
% INPUTS
% required:
%    y       = vector [1 x K] of observed counts
%    lambda  = vector [1 x K] of arrival rate
%    alpha   = detection probability
%    delta   = vector [1 x K] of death probability
%    i       = index to start the summation
%    a, b, f = parameters of initial factor over variable i (the alpha
%              message), which has the form f(s)exp(as + b) for a bounded
%              degree polynomial f
%
% OUTPUTS
%    likelihood = the likelihood, computed as the normalization constant of
%                 the unnormalized marginal of n_i
%
%    a, b, f    = parameters of generating function for the unnormalized
%                 marginal of, which has the form f(s) exp(as + b)
%
%       - a, b are scalars
%       - f is a vector of coefficients, smallest degree first


% MESSAGES AND REPRESENTATION
%
% Messages are defined as:
%
%  gamma(n_i, n_j) := alpha(n_i) * p(n_j, y_{i+1:j} | n_i) --> F_{ij}(s,t)
%
% Then note that we can recover the unnormalized marginal over n_i as:
%
%   p(n_i) = alpha(n_i) * p(y_{i+1:K} | n_i)
%          = \sum_{n_j} gamma(n_i, n_j) 
%
% The joint PGF G_{ij}(s,t) for the gamma messages is defined as:
%
%  G_{ij}(s,t) := \sum_{n_i, n_j} gamma(n_i, n_j) s^{n_i} t^{n_j}
%
% It has the form
%
%  G_{ij}(s,t) = f(s,t) exp(ast + bs + ct + d)
%
% where f is a bounded-degree polynomial in s and t.

K = length(lambda);

if i == K 
   % In this case, there are no additional variables to eliminate. 
   % The unnormalized marginal is equal to the forward message.

   f = f(:); % always return a column vector, to be consistent with other return cases
   loglikelihood = log(sum(f)) + a + b;
   return; 
end

if i < 1 || i > K
    error('Requested index is out of bounds');
end

% Intialize joint PGF over n_i and survivors to next time step 
% from the input PGF over n_i
%
[a, b, c, d, f] = init_survivors(a, b, f, delta(i));

% Now iterate
for k = i+1:K
    [c, d]     = immigrants(c, d, lambda(k));
    [a, c, f]  = observed(y(k), a, c, f, rho);
    
    [d, f] = normalize(d, f);
    
    %last iteration does not have survivors, it just adds immigrants and conditions on observ
    if k <= length(delta)
        [a, b, c, d, f] = survivors(a, b, c, d, f, delta(k));
    end
end

% Finally: set t = 1 and return function f(s) exp(as + b)
a = a + b;     % terms with s
b = c + d;     % terms without s

f = sum(f,2);  % t = 1 --> sum over rows of f  

% Likelihood: evaluate at s=1
loglikelihood = log(sum(f)) + a + b;

end


function [aPrime, bPrime, cPrime, dPrime, fPrime] = init_survivors(a, b, f, delta)
% INIT: starting with initial univariate PGF for alpha message, 
%       produce initial joint PGF using binomial thinning w/ prob. delta 
%
% Input: parameters of initial PGF A(s) = f(s) exp(as + b)
% 
% Output: parameters f, a, b, c, d of joint PGF 
%
%       G(s,t) := A( s * (delta * t + 1 - delta) )
%      
%               = f(s, t) * exp(a*s*t + b*s + c*t + d)  <-- representation
%

aPrime = a *   delta;
bPrime = a * (1 - delta);
cPrime = 0;
dPrime = b;

% Compute polynomial coefficients
%
%    f(s,t) = f(s * (delta * t + 1 - delta)
%
%           = sum_i f(i) s^i (delta * t + 1 - delta)^i
%
%  The ith term in the sum above is a polynomial in t times a fixed power
%  s^i --> this means that the ith row of f(s,t) is equal to f(i) times the
%  coefficients of (delta * t + 1 - delta)^i

degree = length(f) - 1;
fPrime = zeros(degree+1, degree+1);
tPolyBase = [1-delta, delta]; % polynomial (delta*t + 1-delta)
tPoly     = 1;                % ith power of (delta * t + 1 - delta)
for i=1:degree+1
    row = f(i) * tPoly;
    fPrime(i,1:length(row)) = row;
    tPoly = conv(tPoly, tPolyBase);
end
end


function [dPrime, fPrime] = normalize(d, f)

const = max(f(:));
fPrime = f / const;
dPrime = d + log(const);

end


function [cPrime, dPrime] = immigrants(c, d, lambda_k)
% IMMIGRANTS := incorporate the GF of immigrants between t_{k-1} and t_k
%
%  F(s,t) <-- F(s,t)*exp(lambda_k * (t - 1))

cPrime = c + lambda_k;
dPrime = d - lambda_k;

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

function [aPrime, bPrime, cPrime, dPrime, fPrime] = survivors(a, b, c, d, f, delta)
% SURVIVORS := thin the population by delta_k

aPrime = a * delta;
bPrime = b + a * (1 - delta);
cPrime = c * delta;
dPrime = d + c * (1 - delta);

% Perform the composition f'(s,t) = f(s,g(t))
%
%    g(t) = delta_k t + 1 - delta_k
%
% We can do this separately for each row of fPrime, which corresponds to a
% polynomial in t multiplied by a fixed power of s. For each row, use
% Horner's method. Because g(t) is linear, the degree does not increase, so
% we don't need to worry about the dimensions of fPrime increasing.

fPrime = f;
for i=1:size(fPrime,1)
    fPrime(i,:) = compose_poly_horner_special(f(i,:), [1-delta, delta]);
end
end


function [aPrime, cPrime, fPrime] = observed(y, a, c, f, rho)
% OBSERVED := observe y, a value thinned by rho
%
%  New form of PGF is
%             1                    d^y        |
%  F'(s,t) = --- * (t * rho)^y *  ---- F(s,u) |
%             y!                  du^y        | u = t*(1-rho)
%
% We split this into three parts below
%
%  PART 1: Take the yth derivative of F(s,u) / y!
%  PART 2: Evaluate the result at t*(1-rho)
%  PART 3: Multiply the result by (t * rho)^y


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: Take the yth derivative of F(s,u) / y!
%
% The formula for this is
%
%   d^y
%  ---- F(s,u) = exp(asu + bs + cu + d) x
%  du^y
%                          1                        d^x
%          \sum_{x=0}^y --------- (as + c)^{y-x} * ---- f(s,u)
%                       x! (y-x)!                  du^x
%
%          |-------------------- f'(s,u) ---------------------|
%
%
% Note that the exponential part does not change. So we simply need to
% iterate over x and construct the new polynomial part f'.
%
% Other notes:
%  -- The degree wrt s increases by y in this step
%  -- The degree wrt u does not change
%  -- The differentation happens along rows of f, which are polynomials in
%     u (multiplied by a fixed power of s)
%  -- For each value of x, we need to multiply each column of f by the
%     polynomial (as + c)^{y-x}


[m,n] = size(f);

sDegree = m+y;  % s degree increases by y
tDegree = n;

fPrime = zeros(sDegree, tDegree);

xthDerivative = f;
for x = 0:min(y, n-1)
    if x > 0
        %compute the next partial of xthDerivative
        xthDerivative = xthDerivative(:, 2:end); %throw away the constant
        powers = 1:size(xthDerivative,2);
        xthDerivative = bsxfun(@times, xthDerivative, powers);
    end
    
    % Compute (as + c)^{y-x}. This is equal to g(h(s)) for
    %     g(z) := z^{y-x}
    %     h(s) := as + c
    % So, compute it as a polynomial composed with a linear function.
    g = [zeros(1, y-x) 1]; % z^{y-x}
    h = [c a];             % as + c
    sPoly = compose_poly_horner_special(g, h);
    
    % Initialize the bivariate polynomial for this term of the sum
    tDegreeTerm = size(xthDerivative, 2);
    term        = zeros(sDegree, tDegreeTerm);
    
    % Multiply each column of xthDerivative by sPoly
    for j=1:tDegreeTerm
        col                    = conv(xthDerivative(:,j), sPoly);
        term(1:length(col), j) = col;
    end
    
    %compute denominator
    denom = gamma(x+1)*gamma(y-x+1);  % x! * (y-x)!
    
    % Add to running sum
    fPrime = left_aligned_sum(fPrime, term ./ denom);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: Evaluate the result at t*(1-rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform substitution in exponential function
aPrime = a * (1 - rho);
cPrime = c * (1 - rho);

% Perform substitution in polynomial: multiply column j of fPrime by (1-rho)^{j-1}
oneMinusRhoPowers = (1 - rho) .^ ( 0 : tDegree - 1 );
fPrime = bsxfun(@times, fPrime, oneMinusRhoPowers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3: Multiply the result by (t * rho)^y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fPrime = fPrime .* rho^y;             % multiply by rho^y
fPrime = [zeros(sDegree, y) fPrime];  % multiply by t^y

end



function [C] = left_aligned_sum(A,B)
% LEFTALIGNEDSUM := add two vectors of different lengths by padding the shorter vector with zeros

nRows   = size(A,1);

aWidth = size(A,2);
bWidth = size(B,2);

if aWidth < bWidth
    C = [A, zeros(nRows,bWidth - aWidth)] + B;
elseif bWidth < aWidth
    C = [B, zeros(nRows,aWidth - bWidth)] + A;
else
    C = A + B;
end

end