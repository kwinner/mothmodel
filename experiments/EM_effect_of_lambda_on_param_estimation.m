function results = EM_effect_of_lambda_on_param_estimation()

%controllable parameters
lambda_space = [1, 2, 3, 4, 5, 7.5, 10, 12.5, 15, 20, 30];
nRepeats     = 3;

%computed from controlled params
nLambda = numel(lambda_space);

%generally fixed parameters
theta = struct('mu',[],'sigma',[],'lambda',[]);
for i = 1:nLambda
	theta(i).mu     = 4;
	theta(i).sigma  = 2;
	theta(i).lambda = lambda_space(i);
end
theta_cold = struct('mu',1,'sigma',1,'lambda',1);

params = struct('N',500,'alpha',1,'t',1:10);

state = struct('y',[],'q',[],'n',[],'p',[]);

results = struct;

		save('resultsEMLambdaExp.mat', 'results');
		save('resultsEMLambdaExp.mat', 'results');

%begin experiments
for iRepeat = 1:nRepeats
	for iLambda = 1:nLambda
		iExperiment = (iLambda - 1) * nRepeats + iRepeat;
		fprintf('Experiment #%d/%d\n', iExperiment, nRepeats*nLambda);
		%sample a state
		[state(iExperiment).y, state(iExperiment).q, state(iExperiment).n, state(iExperiment).p] = sampleState(params.N, theta(iLambda).mu, theta(iLambda).sigma, theta(iLambda).lambda, params.alpha, params.t);

		%do the different experiments
		%save the params into the result struct for analysis later
		results(iExperiment).iRepeat = iRepeat;
		results(iExperiment).iLambda = iLambda;
		results(iExperiment).theta   = theta(iLambda);
		results(iExperiment).params  = params;
		results(iExperiment).state   = state(iExperiment);

		%y shouldn't be in params, but I don't have time to change it all now...
		%it behaves like a param during EM/mcmc
		params_with_y = params;
		params_with_y.y = state(iExperiment).y;

		%give zonneveld a crack at it
		fprintf('\tZonneveld...\n');
		results(iExperiment).theta_zonn = zonneveldtheta(params_with_y);

		fprintf('\tEM from Zonn...\n');
		results(iExperiment).theta_em_from_zonn = expectationmaximization(results(iExperiment).theta_zonn, params_with_y);

		% fprintf('\tEM from cold...\n');
		% results(iExperiment).theta_em_from_cold = expectationmaximization(theta_cold, params_with_y);

		fprintf('\tEM oracle...\n');
		results(iExperiment).theta_em_oracle    = expectationmaximization_oracle(state(iExperiment), theta_cold, params_with_y);

		save('resultsEMLambdaExp.mat', 'results');
	end
end
end

function theta_struct = theta2struct(theta_array)
theta_struct = struct('mu',theta_array(1),'sigma',theta_array(2),'lambda', theta_array(3));
end