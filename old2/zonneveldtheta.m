function [theta_zonn] = zonneveldtheta(params) 

options = optimoptions('fmincon', 'Algorithm', 'interior-point','Diagnostics','on');%, 'Display', 'on');

problem = struct;
problem.objective = @(theta) -zonneveldLL(theta, params.y, params.N, params.alpha, params.t);
problem.x0        = [1,1,1];
problem.lb        = [1,1,1];
problem.ub        = [50,50,50];
problem.solver    = 'fmincon';
problem.options   = options;

theta_opt = fmincon(problem);
theta_zonn = struct;
theta_zonn.mu = theta_opt(1);
theta_zonn.sigma = theta_opt(2);
theta_zonn.lambda = theta_opt(3);

end