function [theta_zonn] = zonneveldtheta(y, params, varargin)

DEFAULT_THETA_0  = [1 1 1];
DEFAULT_THETA_LB = [1 1 1];
DEFAULT_THETA_UB = [100 100 100];

parser = inputParser;

addParamValue(parser, 'theta_0',  DEFAULT_THETA_0);
addParamValue(parser, 'theta_lb', DEFAULT_THETA_LB);
addParamValue(parser, 'theta_ub', DEFAULT_THETA_UB);

parser.parse(varargin{:});
theta_0  = parser.Results.theta_0;
theta_lb = parser.Results.theta_lb;
theta_ub = parser.Results.theta_ub;

options = optimoptions('fmincon', 'Algorithm', 'interior-point','Display','off');

problem = struct;
problem.objective = @(theta) -zonneveldLL(theta, y, params);
problem.x0        = theta_0;
problem.lb        = theta_lb;
problem.ub        = theta_ub;
problem.solver    = 'fmincon';
problem.options   = options;

theta_opt = fmincon(problem);
theta_zonn = struct;
theta_zonn.mu = theta_opt(1);
theta_zonn.sigma = theta_opt(2);
theta_zonn.lambda = theta_opt(3);

end