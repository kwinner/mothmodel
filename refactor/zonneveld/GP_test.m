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
	problem.objective = @(theta) -GPLL(struct('mu',theta(1),'sigma',theta(2),'lambda',theta(3)),y, struct('N',theta(4),'t',t,'alpha',alpha));
	% problem.x0 = [mu, sigma, lambda, N];
	problem.x0 = [mu sigma lambda N];
	problem.lb = [-5, 1e-6, 1e-6, 1];
	problem.ub = [30, 10, 31, 2000];
	problem.solver = 'fmincon';
	problem.options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

	[theta_zonn,~,~,~,~,~,theta_hessian] = fmincon(problem);

	fprintf('%d: %.2f %.2f %.2f %d\n', iter, theta_zonn(1), theta_zonn(2), theta_zonn(3), theta_zonn(4));

	theta_cov = inv(theta_hessian);
	CI_width = abs(2.*sqrt(diag(theta_cov)))';

	theta_error = abs(theta_zonn - [mu, sigma, lambda, N]);
	coverage = coverage + (theta_error <= CI_width);
end