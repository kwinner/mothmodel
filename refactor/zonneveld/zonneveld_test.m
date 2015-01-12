%zonneveld_test

%parameters learned for population A1 in the Gross et al paper
mu     = 6.9;
sigma  = 3.4;
lambda = 1/0.645;
alpha  = .32;
N      = 216;

realdata = sfs;
t = [realdata(~isnan([realdata.A1])).t] + 1;
y = [realdata(~isnan([realdata.A1])).A1];

coverage = [0,0,0,0];
for iter = 1:100
	y = sampleState(struct('mu',mu,'sigma',sigma,'lambda',lambda), struct('N',N,'t',t,'alpha',alpha));
	problem = struct;
	problem.objective = @(theta) -zonn_loglikelihood(y, theta(1:2), theta(3), t, alpha, theta(4));
	% problem.x0 = [mu, sigma, lambda, N];
	problem.x0 = [1 1 1 1];
	problem.lb = [-inf, 0, 0, 1];
	problem.ub = [inf, inf, inf, inf];
	problem.solver = 'fmincon';
	problem.options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

	[theta_zonn,~,~,~,~,~,theta_hessian] = fmincon(problem);

	fprintf('%.2f %.2f %.2f %d\n', theta_zonn(1), theta_zonn(2), theta_zonn(3), theta_zonn(4));

	theta_cov = inv(theta_hessian);
	CI_width = abs(2.*sqrt(diag(theta_cov)))';

	theta_error = abs(theta_zonn - [mu, sigma, lambda, N]);
	coverage = coverage + (theta_error <= CI_width);
end