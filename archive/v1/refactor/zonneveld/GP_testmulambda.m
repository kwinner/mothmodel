%zonneveld_test

%parameters learned for population A1 in the Gross et al paper
mu     = 9.5;
sigma  = 2.8;
lambda = 7;
alpha  = 1;
N      = 50;

t = 0:5:40;

coverage = [0,0];
theta = cell(100,1);
theta_hessian = cell(100,1);
theta_error = cell(100,1);
CI_width = cell(100,1);
fprintf(['\n' repmat('.',1,100) '\n\n']);
parfor iter = 1:100
	fprintf('\b|\n');
	y = sampleState(struct('mu',mu,'sigma',sigma,'lambda',lambda), struct('N',N,'t',t,'alpha',alpha));
	problem = struct;
	problem.objective = @(theta) -GPLL(struct('mu',theta(1), 'sigma',sigma, 'lambda',theta(2)), y, struct('t', t, 'alpha',alpha, 'N',N));
	% problem.x0 = [mu, sigma, lambda, N];
	problem.x0 = [mu lambda];
	problem.lb = [-5, 1e-6];
	problem.ub = [45, 31];
	problem.solver = 'fmincon';
	problem.options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

	[theta_zonn{iter},~,~,~,~,~,theta_hessian{iter}] = fmincon(problem);

% 	fprintf('%.2f %.2f %.2f %d\n', theta_zonn(1), theta_zonn(2), theta_zonn(3), theta_zonn(4));

	theta_cov = inv(theta_hessian{iter});
	CI_width{iter} = abs(2.*sqrt(diag(theta_cov)))';

	theta_error{iter} = abs(theta_zonn{iter} - [mu, lambda]);
	coverage = coverage + (theta_error{iter} <= CI_width{iter});
end