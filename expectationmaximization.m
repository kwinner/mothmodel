function [mu, sigma, lambda] = expectationmaximization (y, mu, sigma, lambda, N, alpha, T)

N_ITERATIONS = 100;

%options for fminunc
options = optimoptions('fmincon', 'DerivativeCheck', 'on', 'GradObj', 'on');

%the start state for the E step is always the same
[q_0, n_0] = naiveStateExplanation(y,N);

%start the loop
for iter = 1:N_ITERATIONS
	%compute the p matrix
	K = numel(T);
	P = zeros(K+1,K+1);
	for i = 1:K+1
		for j = i:K+1
			P(i,j) = pdf_p(i, j, mu, sigma, lambda, @(i) Tlookup(T,i));
		end
	end

	Q = mcmc(P, q_0, n_0, y, alpha);
	q = Q{end}; %this should be an average

	problem = struct;
	problem.objective = @(theta) objective(theta(1), theta(2), theta(3), @(i)Tlookup(T,i), q);
	problem.x0        = [mu, sigma, lambda];
	problem.lb        = [1, 1, 1];
	problem.solver    = 'fmincon';
	problem.options   = options;
	theta = fmincon(problem);
	mu = theta(1); sigma = abs(theta(2)); lambda = max(theta(3),0);
end

end

function [y, g] =  objective (mu, sigma, lambda, T, Q)
%compute P from mu, sigma, lambda
P = zeros(size(Q));
for i = 1:size(P,1)
	for j = i:size(P,2)
		P(i,j) = pdf_p(i, j, mu, sigma, lambda, T);
	end
end

%the objective function is simply the NLL of mu, sigma, lambda
y = NLLtheta(Q, P);

g = zeros(1,3);

%compute the gradient w.r.t. mu
for i = 1:size(P,1)
	for j = i:size(P,2)
		g(1) = g(1) + gradientWRTmu(i, j, P(i,j), Q(i,j), mu, sigma, lambda, T);
	end
end
%compute the gradient w.r.t. sigma
for i = 1:size(P,1)
	for j = i:size(P,2)
		g(2) = g(2) + gradientWRTsigma(i, j, P(i,j), Q(i,j), mu, sigma, lambda, T);
	end
end
%compute the gradient w.r.t. lambda
for i = 1:size(P,1)
	for j = i:size(P,2)
		g(3) = g(3) + gradientWRTlambda(i, j, P(i,j), Q(i,j), mu, sigma, lambda, T);
	end
end

end

%% utility to put -inf and inf bounds on T while maintaining a sensible indexing scheme
function [t] = Tlookup (T, i)
	if i <= 0
		t = -inf;
	elseif i > numel(T)
		t = inf;
	else
		t = T(i);
	end
end
		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            neg log likelihood of theta given q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%note: for optimization, P is precomputed as a function of theta
function [NLL] = NLLtheta (Q, P)
N = sum(Q(:));

NLL = logfactorial(N);

cellLL = logfactorial(Q) + Q .* logfactorial(P);

%nans occur where P == 0, but correspondingly for those cells Q == 0, so 0 log 0 = 0
NLL = NLL + sum(cellLL(~isnan(cellLL(:))));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           PDF for cell prob as a func of theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_ij] = pdf_p (i, j, mu, sigma, lambda, T)
integrand = @(s) normpdf(s,mu,sigma) .* (expcdf(T(j)-s,lambda) - expcdf(T(j-1)-s,lambda));
p_ij = quadgk(integrand, T(i-1), T(i));
end
% function [p_ij] = pdf_p (i, j, mu, sigma, lambda, T)
% %the part of f_p inside the integral, integrated from t_i to t_i+1
% f_S = @(s) normpdf(s, mu, sigma) .* exp(lambda .* s);

% p_ij = (exp(-lambda * T(j-1)) - exp(-lambda * T(j))) * quadgk(f_S, T(i-1), T(i));

% if p_ij == 0
% 	disp hello
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         dp/dmu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dmu sounds like french for 'the emu'
function [grad] = gradientWRTmu (i, j, p_ij, q_ij, mu, sigma, lambda, T)
%the form of the gradient is always q_i/p_i * dp_i/dtheta_k
grad = q_ij/p_ij;;

%inner gradient: dp_i/dmu
%lambda component:
grad = grad * (exp(-lambda * T(j-1)) - exp(-lambda * T(j)));
%'inner' component:
f_S = @(s) (s-mu)./(sigma^2) .* exp(lambda .* s) .* normpdf(s, mu, sigma);
grad = grad * quadgk(f_S, T(i-1), T(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      dp/dsigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grad] = gradientWRTsigma (i, j, p_ij, q_ij, mu, sigma, lambda, T)
%the form of the gradient is always q_i/p_i * dp_i/dtheta_k
grad = q_ij/p_ij;;

%inner gradient: dp_i/dmu
%lambda component:
grad = grad * (exp(-lambda * T(j-1)) - exp(-lambda * T(j)));
%'inner' component:
f_S = @(s) (((s-mu).^2)./(sigma^3) - 1/sigma) .* exp(lambda .* s) .* normpdf(s, mu, sigma);
grad = grad * quadgk(f_S, T(i-1), T(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     dp/dlambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grad] = gradientWRTlambda (i, j, p_ij, q_ij, mu, sigma, lambda, T)
%the form of the gradient is always q_i/p_i * dp_i/dtheta_k
grad = q_ij/p_ij;;

%inner gradient: dp_i/dmu
%note, the lambda gradient is unfortunately factored into two pieces (first, second)
%first lambda component:
first = (T(j) * exp(-lambda * T(j)) - T(j) * exp(-lambda * T(j-1)));
%first 'inner' component:
f1_S = @(s) exp(lambda .* s) .* normpdf(s, mu, sigma);
first = first * quadgk(f1_S, T(i-1), T(i));

%second lambda component:
second = (exp(-lambda * T(j-1)) - exp(-lambda * T(j)));
%second 'inner' component:
f2_S = @(s) s .* exp(lambda .* s) .* normpdf(s, mu, sigma);
second = second * quadgk(f2_S, T(i-1), T(i));

grad = grad * (first + second);
end