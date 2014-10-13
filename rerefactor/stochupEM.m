function [theta, state, runtime] = stochupEM(y, params, theta_0, varargin)

DEFAULT_EM_ITERATIONS   = 200;           %how many iterations to run EM
DEFAULT_MIXING_PARAM    = .2;            %what percent of Q from E step to add to the running Q
DEFAULT_MCMC_ITERATIONS = 100;           %how many iterations to run MCMC for during E step
DEFAULT_MCMC_BURNIN     = 50;            %how many iterations of MCMC to throw away
DEFAULT_GRADIENTS       = 'off';         %whether to use my gradient functions
DEFAULT_THETA_LB        = [1 1 1];       %lower bounds on the values of theta
DEFAULT_THETA_UB        = [100 100 100]; %upper bounds on the values of theta

parser = inputParser;

addRequired(parser, 'y');
addRequired(parser, 'params');
addRequired(parser, 'theta_0');
addParamValue(parser, 'state_0',         []);
addParamValue(parser, 'EM_iterations',   DEFAULT_EM_ITERATIONS);
addParamValue(parser, 'mixing_param' ,   DEFAULT_MIXING_PARAM);
addParamValue(parser, 'MCMC_iterations', DEFAULT_MCMC_ITERATIONS);
addParamValue(parser, 'MCMC_burnin',     DEFAULT_MCMC_BURNIN);
addParamValue(parser, 'gradients',       DEFAULT_GRADIENTS);
addParamValue(parser, 'theta_lb',        DEFAULT_THETA_LB);
addParamValue(parser, 'theta_ub',        DEFAULT_THETA_UB);
addParamValue(parser, 'oracle',          false);

parser.parse(y, params, theta_0, varargin{:});
y               = parser.Results.y;
params          = parser.Results.params;
theta_0         = parser.Results.theta_0;
state_0         = parser.Results.state_0;
EM_iterations   = parser.Results.EM_iterations;
mixing_param    = parser.Results.mixing_param;
MCMC_iterations = parser.Results.MCMC_iterations;
MCMC_burnin     = parser.Results.MCMC_burnin;
gradients       = parser.Results.gradients;
theta_lb        = parser.Results.theta_lb;
theta_ub        = parser.Results.theta_ub;
oracle          = parser.Results.oracle;

%we do /need/ an initial state, but we can use the naive one if one is not provided
if isempty(state_0)
	state_0 = naiveStateExplanation(y, params.N);
end
if ~isfield(state_0,'p') || isempty(state_0.p)
	state_0.p = ppdf(theta_0, params);
end

%initialize the optimizer options struct (same for all iterations)
if strcmp(gradients, 'on')
	options = optimoptions('fmincon', 'Display', 'notify', 'GradObj', 'on', 'DerivativeCheck', 'off');
else
	options = optimoptions('fmincon', 'Display', 'notify', 'GradObj', 'off', 'Algorithm', 'interior-point');
end

%initialize theta and state, the two arrays that track the progress of the EM
theta = theta_0;
state{1} = state_0;
runtime = 0;
qhat = state_0.q;

%start the main loop
tic
for iter = 2:EM_iterations+1
	fprintf('iteration #%d/%d... ', iter-1, EM_iterations);

	%% ---E STEP---
	if oracle
		%oracle means state_0 should have been the true state, don't need any more E step
		state{iter} = state{iter-1};
	else
		state{iter} = mcmc(y, state{iter-1}(end), params, 'iterations', MCMC_iterations);
	end

	%% -processing-
	%compute the average q from the final runs of the E step
	if length(state{iter}) > 1
		qmean = mean(cat(3,state{iter}(max(1,end-MCMC_burnin):end).q),3);
		qhat = (1-mixing_param) .* qhat + mixing_param .* qmean;
	else
		qhat = (1-mixing_param) .* qhat + mixing_param .* state{iter}.q;
    end

	%% ---M STEP---
	problem = struct;
	problem.objective = @(theta) objective(theta, qhat, params);
	problem.x0        = [theta(iter-1).mu, theta(iter-1).sigma, theta(iter-1).lambda];
	problem.lb        = theta_lb;
	problem.ub        = theta_ub;
	problem.solver    = 'fmincon';
	problem.options   = options;
	thetaArray         = fmincon(problem);
	theta(iter).mu     = thetaArray(1);
	theta(iter).sigma  = thetaArray(2);
	theta(iter).lambda = thetaArray(3);

	%% ------------

	%compute and display some runtime stats
	runtime(iter) = toc - sum(runtime);
	totalruntime = sum(runtime);
	avgruntime = totalruntime / (iter-1);
	estremainingtime = avgruntime * (EM_iterations - iter + 1);
	fprintf('time: %.1f/%.1f, remaining: %.1f\n', runtime(iter), totalruntime, estremainingtime);
end

end

function [y] = objective (theta, qhat, params)
%compute the pdf of p subject to the new values for theta
p_theta = ppdf(theta, params);

y = -sum(qhat(p_theta ~= 0) .* log(p_theta(p_theta ~= 0)));

end