% function [theta, CI] = MLE_theta( y, T, learningmask, mu, sigma, lambda, N, alpha, varargin )
function [theta_hat, CI] = MLE_theta( y, T, varargin );

%start timing
tic

DEFAULT_OBJECTIVE     = 'zonn';
%                              mu         sigma      lambda          N      alpha
DEFAULT_LEARNINGMASK  = [       1,          1,          1,           1,       0  ];
DEFAULT_THETA_LB      = [ min(T)-var(T),  1e-6,       1e-6,          1,      1e-6];
DEFAULT_THETA_UB      = [ max(T)+var(T), var(T),  max(T)-min(T), 100*sum(y),  1  ];
DEFAULT_THETA_0       = [    mean(T),    var(T),     var(T),       sum(y),   0.5 ];

DEFAULT_OPTIM_ALG = 'interior-point';
DEFAULT_DISPLAY   = false;

parser = inputParser;

addOptional  (parser, 'learningmask', DEFAULT_LEARNINGMASK);
addOptional  (parser, 'mu',           []);
addOptional  (parser, 'sigma',        []);
addOptional  (parser, 'lambda',       []);
addOptional  (parser, 'N',            []);
addOptional  (parser, 'alpha',        []);
addParamValue(parser, 'objective',    DEFAULT_OBJECTIVE);
addParamValue(parser, 'optim_alg',    DEFAULT_OPTIM_ALG);
addParamValue(parser, 'display',      DEFAULT_DISPLAY);

parser.parse(varargin{:});
learningmask = parser.Results.learningmask;
mu           = parser.Results.mu;
sigma        = parser.Results.sigma;
lambda       = parser.Results.lambda;
N            = parser.Results.N;
alpha        = parser.Results.alpha;
objective    = parser.Results.objective;
optim_alg    = parser.Results.optim_alg;
display      = parser.Results.display;

options = optimoptions('fmincon', 'Algorithm', optim_alg, 'Display', 'off');

%initialize everything about theta
%note this is a very flexible parameter spec with very little checking, use carefully
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

%construct the problem
problem = struct;
switch objective
	case {'GP', 'gaussian'}
		objectivefun = @gaussian_NLL;
	case {'zonn', 'zonneveld'}
		objectivefun = @zonn_NLL;
end

problem.objective = @(theta) objective_wrapper(y, T, learningmask, theta, fixed_params, objectivefun);
problem.x0        = theta_0;
problem.lb        = theta_lb;
problem.ub        = theta_ub;
problem.solver    = 'fmincon';
problem.options   = options;

%run the optimizer
[theta_hat, ~, ~, ~, ~, ~, theta_hessian] = fmincon(problem);

%compute the CIs
theta_cov = inv(theta_hessian);
CI_width = abs(2.*sqrt(diag(theta_cov)))';
CI = [theta_hat - CI_width; theta_hat + CI_width];

runtime = toc;

if display
	msg = 'theta_hat = [';
	for itheta = 1:numel(theta_hat)
		msg = sprintf('%s%.3f %c %.3f', msg, theta_hat(itheta), 177, CI_width(itheta));
		if itheta < numel(theta_hat)
			msg = [msg, ', '];
		end
	end
	msg = [msg, ']'];
	msg = sprintf('%s finished in %.2fs\n', msg, runtime);

	fprintf(msg);
end

end


function objective = objective_wrapper(y, T, learningmask, theta, fixed_params, objectivefun)

params(logical(learningmask))   = theta;
params(logical(1-learningmask)) = fixed_params;

objective = objectivefun(y, T, params(1), params(2), params(3), params(4), params(5));

end
