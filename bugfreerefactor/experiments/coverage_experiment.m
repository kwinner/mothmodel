function coverage_experiment( T_min, T_max, T_freq, mu, sigma, lambda, N, alpha, niters, out_filename )

%process the inputs if they aren't numeric (generally if this run from cmdline)
if isstr(T_min);  T_min  = str2num(T_min);  end
if isstr(T_max);  T_max  = str2num(T_max);  end
if isstr(T_freq); T_freq = str2num(T_freq); end
if isstr(mu);     mu     = str2num(mu);     end
if isstr(sigma);  sigma  = str2num(sigma);  end
if isstr(lambda); lambda = str2num(lambda); end
if isstr(N);      N      = str2num(N);      end
if isstr(alpha);  alpha  = str2num(alpha);  end
if isstr(niters); niters = str2num(niters); end

verbose = false;

SING_WARNING_STATE = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

T = T_min:T_freq:T_max;

if exist('out_filename','var')
	out_file = fopen(out_filename);
else
	%file id 1 is standard out (so intuitive!)
	out_file = 1;
end

%initialize data storage
theta = [mu, sigma, lambda, N];
y = zeros(niters, numel(T));
theta_zonn    = zeros(niters, 4);
CI_width_zonn = zeros(niters, 4);
runtime_zonn  = zeros(niters, 1);
coverage_zonn = zeros(niters, 4);
error_zonn    = zeros(niters, 4);
theta_gaussian    = zeros(niters, 4);
CI_width_gaussian = zeros(niters, 4);
runtime_gaussian  = zeros(niters, 1);
coverage_gaussian = zeros(niters, 4);
error_gaussian    = zeros(niters, 4);

for iter = 1:niters
	%generate a sample for this iteration
	y(iter,:) = samplestate(mu, sigma, lambda, N, alpha, T);

	%learn MLE estimates w/ Zonneveld
	tic
	[theta_zonn(iter,:), CI_zonn]         = MLE_theta(y(iter,:), T, [1 1 1 1 0], mu, sigma, lambda, N, alpha, 'objective', 'zonn',     'display', false);
	runtime_zonn(iter) = toc;

	%annnd repeat with GP
	tic
	[theta_gaussian(iter,:), CI_gaussian] = MLE_theta(y(iter,:), T, [1 1 1 1 0], mu, sigma, lambda, N, alpha, 'objective', 'gaussian', 'display', false);
	runtime_gaussian(iter) = toc;

	%compute coverage, error
	CI_width_zonn(iter,:)     = abs(theta_zonn(iter,:) - CI_zonn(1,:));
	coverage_zonn(iter,:)     = theta >= CI_zonn(1,:)     & theta <= CI_zonn(2,:);
	error_zonn(iter,:)        = abs(theta - theta_zonn(iter,:));
	CI_width_gaussian(iter,:) = abs(theta_gaussian(iter,:) - CI_gaussian(1,:));
	coverage_gaussian(iter,:) = theta >= CI_gaussian(1,:) & theta <= CI_gaussian(2,:);
	error_gaussian(iter,:)    = abs(theta - theta_gaussian(iter,:));

	if verbose
		format long g
		fprintf(out_file, 'iteration %d:\n', iter)
		fprintf(out_file, 'theta = [%.2f, %.2f, %.2f, %.1f]\n', mu, sigma, lambda, N)
		fprintf(out_file, 'y = %s\n', mat2str(y(iter,:)))

		fprintf(out_file, 'zonneveld:\n')
		fprintf(out_file, 'theta_zonn = [%.2f %c %.2f, %.2f %c %.2f, %.2f %c %.2f, %.1f %c %.1f]\n', theta_zonn(iter,1), 177, CI_width_zonn(iter,1), ...
		                                                                                             theta_zonn(iter,2), 177, CI_width_zonn(iter,2), ...
		                                                                                             theta_zonn(iter,3), 177, CI_width_zonn(iter,3), ...
		                                                                                             theta_zonn(iter,4), 177, CI_width_zonn(iter,4))
		fprintf(out_file, 'coverage_zonn = [%d %d %d %d]\n', coverage_zonn(iter,1), coverage_zonn(iter,2), coverage_zonn(iter,3), coverage_zonn(iter,4))
		fprintf(out_file, 'error_zonn = [%.2f, %.2f, %.2f, %.1f]\n', error_zonn(iter,1), error_zonn(iter,2), error_zonn(iter,3), error_zonn(iter,4))
		fprintf(out_file, 'runtime_zonn = %.2fs\n', runtime_zonn(iter))

		fprintf(out_file, 'gaussian:\n')
		fprintf(out_file, 'theta_gaussian = [%.2f %c %.2f, %.2f %c %.2f, %.2f %c %.2f, %.1f %c %.1f]\n', theta_gaussian(iter,1), 177, CI_width_gaussian(iter,1), ...
		                                                                                                 theta_gaussian(iter,2), 177, CI_width_gaussian(iter,2), ...
		                                                                                                 theta_gaussian(iter,3), 177, CI_width_gaussian(iter,3), ...
		                                                                                                 theta_gaussian(iter,4), 177, CI_width_gaussian(iter,3))
		fprintf(out_file, 'coverage_gaussian = [%d %d %d %d]\n', coverage_gaussian(iter,1), coverage_gaussian(iter,2), coverage_gaussian(iter,3), coverage_gaussian(iter,4))
		fprintf(out_file, 'error_gaussian = [%.2f, %.2f, %.2f, %.1f]\n', error_gaussian(iter,1), error_gaussian(iter,2), error_gaussian(iter,3), error_gaussian(iter,4))
		fprintf(out_file, 'runtime_gaussian = %.2fs\n', runtime_gaussian(iter))

		fprintf(out_file, '\n\n')
	end
end

%compute summary statistics
mean_theta_zonn        = mean(theta_zonn, 1);
mean_CI_width_zonn     = mean(CI_width_zonn, 1);
mean_error_zonn        = mean(error_zonn, 1);
mean_coverage_zonn     = mean(coverage_zonn, 1);
mean_runtime_zonn      = mean(runtime_zonn);

mean_theta_gaussian    = mean(theta_gaussian, 1);
mean_CI_width_gaussian = mean(CI_width_gaussian, 1);
mean_error_gaussian    = mean(error_gaussian, 1);
mean_coverage_gaussian = mean(coverage_gaussian, 1);
mean_runtime_gaussian  = mean(runtime_gaussian, 1);

%print results
fprintf(out_file, 'zonneveld:\n')
fprintf(out_file, 'mean_theta_zonn = [%.2f, %.2f, %.2f, %.1f]\n', mean_theta_zonn(1), mean_theta_zonn(2), mean_theta_zonn(3), mean_theta_zonn(4))
fprintf(out_file, 'mean_CI_width_zonn = [%.2f, %.2f, %.2f, %.1f]\n', mean_CI_width_zonn(1), mean_CI_width_zonn(2), mean_CI_width_zonn(3), mean_CI_width_zonn(4))
fprintf(out_file, 'mean_error_zonn = [%.2f, %.2f, %.2f, %.1f]\n', mean_error_zonn(1), mean_error_zonn(2), mean_error_zonn(3), mean_error_zonn(4))
fprintf(out_file, 'mean_coverage_zonn = [%.2f, %.2f, %.2f, %.2f]\n', mean_coverage_zonn(1), mean_coverage_zonn(2), mean_coverage_zonn(3), mean_coverage_zonn(4))
fprintf(out_file, 'mean_runtime_zonn = %.2fs\n', mean_runtime_zonn)
fprintf(out_file, '\n')
fprintf(out_file, 'gaussian:\n')
fprintf(out_file, 'mean_theta_gaussian = [%.2f, %.2f, %.2f, %.1f]\n', mean_theta_gaussian(1), mean_theta_gaussian(2), mean_theta_gaussian(3), mean_theta_gaussian(4))
fprintf(out_file, 'mean_CI_width_gaussian = [%.2f, %.2f, %.2f, %.1f]\n', mean_CI_width_gaussian(1), mean_CI_width_gaussian(2), mean_CI_width_gaussian(3), mean_CI_width_gaussian(4))
fprintf(out_file, 'mean_error_gaussian = [%.2f, %.2f, %.2f, %.1f]\n', mean_error_gaussian(1), mean_error_gaussian(2), mean_error_gaussian(3), mean_error_gaussian(4))
fprintf(out_file, 'mean_coverage_gaussian = [%.2f, %.2f, %.2f, %.2f]\n', mean_coverage_gaussian(1), mean_coverage_gaussian(2), mean_coverage_gaussian(3), mean_coverage_gaussian(4))
fprintf(out_file, 'mean_runtime_gaussian = %.2fs\n', mean_runtime_gaussian)

warning(SING_WARNING_STATE.state, 'MATLAB:nearlySingularMatrix');

end

