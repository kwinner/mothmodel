function [theta_opt] = zonneveldtheta(params) 

options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

problem = struct;
problem.objective = @(theta) -zonneveldLL(theta, params.y, params.N, params.alpha, params.t);
problem.x0        = [1,1,1];
problem.lb        = [1,1,1];
problem.solver    = 'fmincon';
problem.options   = options;

theta_opt = fmincon(problem);

end