function coverage = CIcoverage(theta_0, params, varargin)

DEFAULT_ITERATIONS = 2;

parser = inputParser;

addRequired  (parser, 'theta_0');
addRequired  (parser, 'params');
addParamValue(parser, 'iterations', DEFAULT_ITERATIONS);

parser.parse(theta_0, params, varargin{:});

theta_0     = parser.Results.theta_0;
params      = parser.Results.params;
niterations = parser.Results.iterations;

theta_0_mat = [theta_0.mu, theta_0.sigma, theta_0.lambda];

nparams   = numel(theta_0_mat);
coverage  = zeros(niterations, nparams);

for iteration = 1:niterations
	% if mod(iteration, 10) == 0
	% 	fprintf('iteration #%d...\n', iteration);
	% end
	
	%sample from the true theta
	y_i     = sampleState(theta_0, params);

	%learn theta from sample, including the CI
	[theta_i, CI_i] = zonneveldtheta(y_i, params);

	%check which CIs cover theta_0
	for iparam = 1:nparams
		if CI_i(1, iparam) <= theta_0_mat(iparam) && CI_i(2, iparam) >= theta_0_mat(iparam)
			coverage(iteration, iparam) = 1;
		end
	end
end

%normalize the coverage counter
coverage = sum(coverage, 1);
coverage = coverage ./ niterations;

end