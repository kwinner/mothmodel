function error = errorofresults(results)

error = struct;

true_thetas = [results.theta];
lambda_values = unique([true_thetas.lambda]);
for iLambda = 1:numel(lambda_values)
	lambda_value = lambda_values(iLambda);
	error(iLambda).lambda = lambda_value; %for the x axis
	error(iLambda).theta_zonn = 0;
	error(iLambda).theta_em_from_zonn = 0;
	error(iLambda).theta_em_oracle = 0;
	error(iLambda).theta_zonn_rel = 0;
	error(iLambda).theta_em_from_zonn_rel = 0;
	error(iLambda).theta_em_oracle_rel = 0;
	error(iLambda).theta_zonn_mu = 0;
	error(iLambda).theta_em_from_zonn_mu = 0;
	error(iLambda).theta_em_oracle_mu = 0;
	error(iLambda).theta_zonn_sigma = 0;
	error(iLambda).theta_em_from_zonn_sigma = 0;
	error(iLambda).theta_em_oracle_sigma = 0;
	error(iLambda).theta_zonn_lambda = 0;
	error(iLambda).theta_em_from_zonn_lambda = 0;
	error(iLambda).theta_em_oracle_lambda = 0;
	
	results_lambda = results([true_thetas.lambda] == lambda_value);
	nRepeats = numel(results_lambda);
	for iRepeat = 1:nRepeats
		error(iLambda).theta_zonn         = error(iLambda).theta_zonn         + rmse(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_zonn) / nRepeats;
		error(iLambda).theta_em_from_zonn = error(iLambda).theta_em_from_zonn + rmse(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_em_from_zonn(end)) / nRepeats;
		error(iLambda).theta_em_oracle    = error(iLambda).theta_em_oracle    + rmse(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_em_oracle(end)) / nRepeats;

		error(iLambda).theta_zonn_rel         = error(iLambda).theta_zonn_rel         + relativeerror(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_zonn) / nRepeats;
		error(iLambda).theta_em_from_zonn_rel = error(iLambda).theta_em_from_zonn_rel + relativeerror(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_em_from_zonn(end)) / nRepeats;
		error(iLambda).theta_em_oracle_rel    = error(iLambda).theta_em_oracle_rel    + relativeerror(results_lambda(iRepeat).theta, results_lambda(iRepeat).theta_em_oracle(end)) / nRepeats;

		error(iLambda).theta_zonn_mu         = error(iLambda).theta_zonn_mu         + rmse(results_lambda(iRepeat).theta.mu, results_lambda(iRepeat).theta_zonn.mu) / nRepeats;
		error(iLambda).theta_em_from_zonn_mu = error(iLambda).theta_em_from_zonn_mu + rmse(results_lambda(iRepeat).theta.mu, results_lambda(iRepeat).theta_em_from_zonn(end).mu) / nRepeats;
		error(iLambda).theta_em_oracle_mu    = error(iLambda).theta_em_oracle_mu    + rmse(results_lambda(iRepeat).theta.mu, results_lambda(iRepeat).theta_em_oracle(end).mu) / nRepeats;

		error(iLambda).theta_zonn_sigma         = error(iLambda).theta_zonn_sigma         + rmse(results_lambda(iRepeat).theta.sigma, results_lambda(iRepeat).theta_zonn.sigma) / nRepeats;
		error(iLambda).theta_em_from_zonn_sigma = error(iLambda).theta_em_from_zonn_sigma + rmse(results_lambda(iRepeat).theta.sigma, results_lambda(iRepeat).theta_em_from_zonn(end).sigma) / nRepeats;
		error(iLambda).theta_em_oracle_sigma    = error(iLambda).theta_em_oracle_sigma    + rmse(results_lambda(iRepeat).theta.sigma, results_lambda(iRepeat).theta_em_oracle(end).sigma) / nRepeats;

		error(iLambda).theta_zonn_lambda         = error(iLambda).theta_zonn_lambda         + rmse(results_lambda(iRepeat).theta.lambda, results_lambda(iRepeat).theta_zonn.lambda) / nRepeats;
		error(iLambda).theta_em_from_zonn_lambda = error(iLambda).theta_em_from_zonn_lambda + rmse(results_lambda(iRepeat).theta.lambda, results_lambda(iRepeat).theta_em_from_zonn(end).lambda) / nRepeats;
		error(iLambda).theta_em_oracle_lambda    = error(iLambda).theta_em_oracle_lambda    + rmse(results_lambda(iRepeat).theta.lambda, results_lambda(iRepeat).theta_em_oracle(end).lambda) / nRepeats;
	end
end

end