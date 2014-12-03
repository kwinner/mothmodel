function [theta, CI] = MLtheta(y, params, varargin)

if isfield(params, 'N')
	DEFAULT_THETA_0  = [1, 1, 1];
	DEFAULT_THETA_LB = [1 1 1];
	DEFAULT_THETA_UB = [15 10 10];
else
	DEFAULT_THETA_0  = [1, 1, 1, 200];
	DEFAULT_THETA_LB = [1, 1, 1, 1];
	DEFAULT_THETA_UB = [15,10,10,5000];
end

DEFAULT_OBJECTIVE = 'GP';      %GP, zonn

parser = inputParser;

addParamValue(parser, 'theta_0',   DEFAULT_THETA_0);
addParamValue(parser, 'theta_lb',  DEFAULT_THETA_LB);
addParamValue(parser, 'theta_ub',  DEFAULT_THETA_UB);
addParamValue(parser, 'objective', DEFAULT_OBJECTIVE);

parser.parse(varargin{:});
theta_0  = parser.Results.theta_0;
theta_lb = parser.Results.theta_lb;
theta_ub = parser.Results.theta_ub;
objective = parser.Results.objective;

options = optimoptions('fmincon', 'Algorithm', 'interior-point','Display','off');

problem = struct;
switch objective
	case 'GP'
		objectivefun = @GPLL;
	case 'zonn'
		objectivefun = @zonneveldLL;
end
if isfield(params, 'N')
	problem.objective = @(theta) -objectivefun(theta, y, params);
else
	problem.objective = @(theta) -objectivefun(theta, y, struct('alpha', params.alpha, 't', params.t, 'N', theta(4)));
end
problem.x0        = theta_0;
problem.lb        = theta_lb;
problem.ub        = theta_ub;
problem.solver    = 'fmincon';
problem.options   = options;

[theta_mat,~,~,~,~,~,theta_hessian] = fmincon(problem);
theta = struct;
theta.mu = theta_mat(1);
theta.sigma = theta_mat(2);
theta.lambda = theta_mat(3);

%compute the confidence intervals for the learned parameters
theta_cov = inv(theta_hessian);
CI_width = abs(2.*sqrt(diag(theta_cov)))';
CI = [theta_mat - CI_width; theta_mat + CI_width];

end