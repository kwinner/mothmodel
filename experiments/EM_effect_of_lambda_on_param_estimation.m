function results = EM_effect_of_lambda_on_param_estimation()

%controllable parameters
% lambda_space = [2, 4, 8, 12, 16, 20, 30, 50];
lambda_space = [14];
nRepeats     = 1;

%computed from controlled params
nLambda = numel(lambda_space);

%generally fixed parameters
theta = struct('mu',[],'sigma',[],'lambda',[]);
for i = 1:nLambda
	theta(i).mu     = 30;
	theta(i).sigma  = 10;
	theta(i).lambda = lambda_space(i);
end
theta_cold = struct('mu',1,'sigma',1,'lambda',1);

params = struct('N',25,'alpha',.8,'t',10:10:100);

state = struct('p', [], 'q',[],'n',[]);

results = struct;

%begin experiments
for iRepeat = 1:nRepeats
	for iLambda = 1:nLambda
		iExperiment = (iLambda - 1) * nRepeats + iRepeat;
		fprintf('Experiment #%d/%d\n', iExperiment, nRepeats*nLambda);
		fprintf('theta: mu = %.2f, sigma = %.2f, lambda = %.2f\n', theta(iLambda).mu, theta(iLambda).sigma, theta(iLambda).lambda);
		%sample a state
		[y, state(iExperiment)] = sampleState(theta(iLambda), params);

		%do the different experiments
		%save the params into the result struct for analysis later
		results(iExperiment).iRepeat = iRepeat;
		results(iExperiment).iLambda = iLambda;
		results(iExperiment).theta   = theta(iLambda);
		results(iExperiment).params  = params;
		results(iExperiment).state   = state(iExperiment);

		%give zonneveld a crack at it
		fprintf('\tZonneveld...\n');
		results(iExperiment).theta_zonn = zonneveldtheta(y, params);
		fprintf('\tZonneveld theta: mu = %.2f, sigma = %.2f, lambda = %.2f\n', results(iExperiment).theta_zonn.mu, results(iExperiment).theta_zonn.sigma, results(iExperiment).theta_zonn.lambda);

		fprintf('\tEM oracle...\n');
		results(iExperiment).theta_em_oracle = stochupEM(y, params, theta_cold, 'state_0', state(iExperiment), 'EM_iterations', 1);
		fprintf('\tEM oracle theta: mu = %.2f, sigma = %.2f, lambda = %.2f\n', results(iExperiment).theta_em_oracle(end).mu, results(iExperiment).theta_em_oracle(end).sigma, results(iExperiment).theta_em_oracle(end).lambda);

		fprintf('\tEM from Zonn...\n');
		results(iExperiment).theta_em_from_zonn = stochupEM(y, params, results(iExperiment).theta_zonn, 'EM_iterations', 20);
		fprintf('\tEM from Zonn theta: mu = %.2f, sigma = %.2f, lambda = %.2f\n', results(iExperiment).theta_em_from_zonn(end).mu, results(iExperiment).theta_em_from_zonn(end).sigma, results(iExperiment).theta_em_from_zonn(end).lambda);

		save('resultsEMLambdaExp.mat', 'results');
	end
end
end