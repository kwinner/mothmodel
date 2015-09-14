function [theta_hat_history, EQ_history, runtime_history] = stochupEM( y, T, varargin )

DEFAULT_EM_ITERATIONS   = 200;           %how many iterations to run EM
DEFAULT_MIXING_PARAM    = .2;            %what percent of Q from E step to add to the running Q
DEFAULT_INITIAL_MCMC_ITERATIONS = 50000; %how many iterations to run for the first iteration of MCMC
DEFAULT_MCMC_ITERATIONS = 100;           %how many iterations to run MCMC for during E step
DEFAULT_INITIAL_MCMC_BURNIN = 40000;      %how many iterations of the first MCMC pass to throw away
DEFAULT_MCMC_BURNIN     = 50;            %how many iterations of MCMC to throw away

%                              mu         sigma      lambda          N      alpha
DEFAULT_LEARNINGMASK  = [       1,          1,          1,           1,       0  ];
DEFAULT_THETA_LB      = [ min(T)-var(T),  1e-6,       1e-6,          1,      1e-6];
DEFAULT_THETA_UB      = [ max(T)+var(T), var(T),  max(T)-min(T), 100*sum(y),  1  ];
DEFAULT_THETA_0       = [    mean(T),    var(T),     var(T),       sum(y),   0.5 ];

parser = inputParser;

addOptional  (parser, 'learningmask',            DEFAULT_LEARNINGMASK);
addOptional  (parser, 'mu',                      []);
addOptional  (parser, 'sigma',                   []);
addOptional  (parser, 'lambda',                  []);
addOptional  (parser, 'N',                       []);
addOptional  (parser, 'alpha',                   []);
addParamValue(parser, 'Q_0',                     []);
addParamValue(parser, 'EM_iterations',           DEFAULT_EM_ITERATIONS);
addParamValue(parser, 'mixing_param' ,           DEFAULT_MIXING_PARAM);
addParamValue(parser, 'initial_MCMC_iterations', DEFAULT_INITIAL_MCMC_ITERATIONS);
addParamValue(parser, 'MCMC_iterations',         DEFAULT_MCMC_ITERATIONS);
addParamValue(parser, 'initial_MCMC_burnin',     DEFAULT_INITIAL_MCMC_BURNIN);
addParamValue(parser, 'MCMC_burnin',             DEFAULT_MCMC_BURNIN);
addParamValue(parser, 'oracle',                  false);
addParamValue(parser, 'use_ARS',                 false);

parser.parse(varargin{:});
learningmask            = parser.Results.learningmask;
mu                      = parser.Results.mu;
sigma                   = parser.Results.sigma;
lambda                  = parser.Results.lambda;
N                       = parser.Results.N;
alpha                   = parser.Results.alpha;
Q_0                     = parser.Results.Q_0;
EM_iterations           = parser.Results.EM_iterations;
mixing_param            = parser.Results.mixing_param;
initial_MCMC_iterations = parser.Results.initial_MCMC_iterations;
MCMC_iterations         = parser.Results.MCMC_iterations;
initial_MCMC_burnin     = parser.Results.initial_MCMC_burnin;
MCMC_burnin             = parser.Results.MCMC_burnin;
oracle                  = parser.Results.oracle;
use_ARS                 = parser.Results.use_ARS;

%parse the (obnoxiously complicated) values for theta
parameters = {mu, sigma, lambda, N, alpha};

if numel(learningmask) ~= 5
	learningmask = DEFAULT_LEARNINGMASK;
end

theta_0      = [];
theta_lb     = [];
theta_ub     = [];
fixed_params = [];

for itheta = 1:5
	if learningmask(itheta)
		%parameter i is to be learned
		theta_lb_i = DEFAULT_THETA_LB(itheta);
		theta_ub_i = DEFAULT_THETA_UB(itheta);
		theta_0_i  = DEFAULT_THETA_0(itheta);
		switch numel(parameters{itheta})
		case 0
			%use all defaults
			[]; %NOOP
		case 1
			%override theta_0
			theta_0_i = parameters{itheta};
		case 2
			%override bounds
			theta_lb_i = min(parameters{itheta});
			theta_ub_i = max(parameters{itheta});
		case 3
			%override all defaults
			theta_lb_i = min(parameters{itheta});
			theta_ub_i = max(parameters{itheta});
			theta_0_i  = median(parameters{itheta});
		end

		theta_lb = [theta_lb, theta_lb_i];
		theta_ub = [theta_ub, theta_ub_i];
		theta_0  = [theta_0,  theta_0_i];
	else
		%parameter i is fixed
		if numel(parameters{itheta}) == 0
			fixed_params = [fixed_params, DEFAULT_THETA_0(itheta)];
		else
			fixed_params = [fixed_params, parameters{itheta}];
		end
	end
end

%construct a start state if one did not exist previously
if isempty(Q_0)
	Q_0 = naiveStateExplanation(y, theta_0(4));
end

K = numel(T);

theta_hat_history     = zeros(EM_iterations,sum(learningmask)); theta_hat_history(1,:) = theta_0;
EQ_history            = zeros(K+1,K+1,EM_iterations); EQ_history(:,:,1) = Q_0
runtime_history       = zeros(EM_iterations,1);

options = optimoptions('fmincon', 'Display', 'notify', 'GradObj', 'off', 'Algorithm', 'interior-point');

%start the main loop
for iter = 2:EM_iterations
	starttime = tic;
	%% ---E STEP---
	if oracle
		EQ_history(:,:,iter) = EQ_history(:,:,iter-1);
	else
		if iter == 2 %'first' iteration
			nMCMCIterations = initial_MCMC_iterations;
			nMCMCBurnin     = initial_MCMC_burnin;
		else
			nMCMCIterations = MCMC_iterations;
			nMCMCBurnin     = MCMC_burnin;
		end

		P = birthdeath_pmf(theta_hat_history(iter-1,1), theta_hat_history(iter-1,2), theta_hat_history(iter-1,3), theta_hat_history(iter-1,4), T);

		EQ_history(:,:,iter) = mcmc(y, ...
			                        EQ_history(:,:,iter-1), ...
			                        theta_hat_history(iter-1,4), ...
			                        theta_hat_history(iter-1,5), ...
			                        'nIterations', nMCMCIterations, ...
			                        'useARS', use_ARS, ...
			                        'burnin', nMCMCBurnin);
	end

	%% ---M STEP---
	problem = struct;
	problem.objective = @(theta) objective(y, T, EQ_history(:,:,iter), learningmask, theta_hat_history(iter-1,:), fixed_params);
	problem.x0        = theta_hat_history(iter-1,:);
	problem.lb        = theta_lb;
	problem.ub        = theta_ub;
	problem.solver    = 'fmincon';
	problem.options   = options;

	theta_hat_history(iter,:) = fmincon(problem);
	runtime_history(iter) = toc(starttime);
end

end


function objective = objective(y, T, EQ, learningmask, theta, fixed_params)

params(logical(learningmask))   = theta;
params(logical(1-learningmask)) = fixed_params;

P = birthdeath_pmf(params(1), params(2), params(3), params(4), T);

objective = -latentvar_LL(EQ,P) - latentvar_obs_LL(y, abundance(y), params(5));

end

