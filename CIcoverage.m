function coverage = CIcoverage(theta_0, params, varargin)

DEFAULT_ITERATIONS = 10;

parser = inputParser;

addRequired  (parser, 'theta_0');
addRequired  (parser, 'params');
addParamValue(parser, 'iterations', DEFAULT_ITERATIONS);

parser.parse(theta_0, params, varargin{:});

theta_0     = parser.Results.theta_0;
params      = parser.Results.params;
niterations = parser.Results.iterations;

theta_0_mat = [theta_0.mu, theta_0.sigma, theta_0.lambda];

nparams       = numel(theta_0_mat);
coverage_GP   = zeros(niterations, nparams);
coverage_zonn = zeros(niterations, nparams);

for iteration = 1:niterations
	% if mod(iteration, 10) == 0
		fprintf('iteration #%d...\n', iteration);
	% end

	%sample from the true theta
	y_i     = sampleState(theta_0, params);

	%learn theta from sample, including the CI
	[theta_GP, CI_GP]     = MLtheta(y_i, params, 'theta_0', theta_0_mat, 'objective', 'GP');
	[theta_zonn, CI_zonn] = MLtheta(y_i, params, 'theta_0', theta_0_mat, 'objective', 'zonn');

	%check which CIs cover theta_0
	for iparam = 1:nparams
		if CI_GP(1, iparam) <= theta_0_mat(iparam) && CI_GP(2, iparam) >= theta_0_mat(iparam)
			coverage_GP(iteration, iparam) = 1;
		end
		if CI_zonn(1, iparam) <= theta_0_mat(iparam) && CI_zonn(2, iparam) >= theta_0_mat(iparam)
			coverage_zonn(iteration, iparam) = 1;
		end
	end
end

%normalize the coverage counter
coverage_GP   = sum(coverage_GP, 1);
coverage_GP   = coverage_GP ./ niterations;
coverage_zonn = sum(coverage_zonn, 1);
coverage_zonn = coverage_zonn ./ niterations;

coverage.GP   = coverage_GP;
coverage.zonn = coverage_zonn;

end