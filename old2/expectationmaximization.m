function [theta, state] = expectationmaximization (theta, params)
%state  = the matrices/vectors making up the observations and outcomes
%		  state.q
%		  state.n (note: this is a function of q)
%theta  = the parameters being optimized by EM
%		  mu, sigma, lambda
%params = the constant parameters
%		  y, N, alpha, t

MAX_EM_ITERATIONS = 200;

%options for fminunc
options = optimoptions('fmincon', 'DerivativeCheck', 'off', 'GradObj', 'on', 'Display', 'off');

[naiveQ, naiveN] = naiveStateExplanation(params.y, params.N);

%initialize the state struct
state = struct;
state.q_samples = cell(1,1);
state.q_samples{1} = naiveQ;

%start the loop
for iter = 2:MAX_EM_ITERATIONS+1
	% if mod(iter-1,5) == 0
		fprintf('\t\titeration #%d/%d\n',iter-1,MAX_EM_ITERATIONS);
	% end
	%% E step (sample new state: q and n)
	%compute p from theta for mcmc
	state(iter-1).p = pdf_p([theta(iter-1).mu, theta(iter-1).sigma, theta(iter-1).lambda],state(iter-1).q_samples{end}, params);
	state(iter).q_samples = mcmc(state(iter-1).p, naiveQ, naiveN, params.y, params.alpha, 'iterations', 2000);
	state(iter).q_avg = mean(cat(3,state(iter).q_samples{end-20:end}),3);
    % state(iter).q = q_samples{end};
	% state(iter).n_avg = abundancy(state(iter).q_avg);
	% state(iter) = state(iter-1);

	%% M step (optimize theta: mu, sigma, lambda)
	problem = struct;
	problem.objective = @(theta) objective(theta, state(iter), params);
	problem.x0        = [theta(iter-1).mu, theta(iter-1).sigma, theta(iter-1).lambda];
	problem.lb        = [1, 1, 1];
	problem.solver    = 'fmincon';
	problem.options   = options;
	theta_packed = fmincon(problem);
	theta(iter).mu     = theta_packed(1);
	theta(iter).sigma  = theta_packed(2);
	theta(iter).lambda = theta_packed(3);
end

end

function [y, g] =  objective (theta, state, params)
%compute p given the new values for theta
state.p = pdf_p(theta, state, params);

%the objective function is simply the NLL of mu, sigma, lambda
%which, if you drop the terms that are ind. of theta is just
%the neg sum of all q_i * log p_i
y = -sum(state.q_avg(state.p~=0) .* log(state.p(state.p~=0)));

%compute the gradients
g = zeros(1,3);
g(1) = gradientWRTtheta(theta, state, params, @gradientWRTmu_integrand);
g(2) = gradientWRTtheta(theta, state, params, @gradientWRTsigma_integrand);
g(3) = gradientWRTtheta(theta, state, params, @gradientWRTlambda_integrand);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           PDF for cell prob as a func of theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = pdf_p (theta, state, params)
%for integration, use -inf and inf as min and max vals on t
T = numel(params.t);
params.t = [-inf, params.t, inf];
P = zeros(T+1);
for i = 1:T+1
	for j = i:T+1
		P(i,j) = quadgk(@(s) pdf_p_integrand(i, j, s, theta, params), params.t(i), params.t(i+1));
	end
end
end
function [integrand] = pdf_p_integrand (i, j, s, theta, params)
mu = theta(1); sigma = theta(2); lambda = theta(3);
[zmin, zmax] = lifespan_domain(params.t(j), params.t(j+1), s);
integrand = normpdf(s, mu, sigma) .* (expcdf(zmax, lambda) - expcdf(zmin, lambda));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      dp/dtheta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%all of the gradients use the same common form, except the form of the integrand
function [grad] = gradientWRTtheta (theta, state, params, integrand_fun)
%for integration, use -inf and inf as min and max vals on t
T = numel(params.t);
params.t = [-inf, params.t, inf];

%grad = - sum(q_i/p_i * d/dmu p_i)
grad = 0;
for i = 1:T+1
	for j = i:T+1
		grad = grad - state.q_avg(i,j) / state.p(i,j) * quadgk(@(s) integrand_fun(i, j, s, theta, params), params.t(i), params.t(i+1));
	end
end
end

%% d/dmu p_i
function [integrand] = gradientWRTmu_integrand (i, j, s, theta, params)
mu = theta(1); sigma = theta(2); lambda = theta(3);
[zmin, zmax] = lifespan_domain(params.t(j), params.t(j+1), s);
integrand = (s-mu)./(sigma.^2) .* normpdf(s, mu, sigma) .* (expcdf(zmax, lambda) - expcdf(zmin, lambda));
end

%% d/dsigma p_i
function [integrand] = gradientWRTsigma_integrand (i, j, s, theta, params)
mu = theta(1); sigma = theta(2); lambda = theta(3);
[zmin, zmax] = lifespan_domain(params.t(j), params.t(j+1), s);
integrand = ((s-mu).^2./(sigma.^3) - 1/sigma) .* normpdf(s, mu, sigma) .* (expcdf(zmax, lambda) - expcdf(zmin, lambda));
end

%% d/dlambda p_i
function [integrand] = gradientWRTlambda_integrand (i, j, s, theta, params)
mu = theta(1); sigma = theta(2); lambda = theta(3);
[zmin, zmax] = lifespan_domain(params.t(j), params.t(j+1), s);
integrand = normpdf(s, mu, sigma);
lhs = zmax ./ lambda .* exppdf(zmax, lambda);
lhs(abs(zmax) == inf) = 0;
rhs = zmin ./ lambda .* exppdf(zmin, lambda);
rhs(abs(zmin) == inf) = 0;
integrand = integrand .* (lhs - rhs);
integrand = -integrand;
end

%% utility function
% there's a problem with the integration of F(zmax) - F(zmin)
% F(z) is 0 for z <= 0, even though z \in [0,inf)
% this isn't a problem for the pdf, but the gradients need to restrict z appropriately
% especially dp/dlambda, which includes multiplication by zmin/zmax
function [zmin, zmax] = lifespan_domain (death_min, death_max, birth)
zmin = death_min - birth;
zmin(death_min == birth) = 0;
zmin = max(zmin, 0); %clamp z to 0
zmax = death_max - birth;
zmax(death_max == birth) = inf;
end