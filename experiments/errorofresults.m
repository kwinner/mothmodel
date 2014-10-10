function error = errorofresults(results)

error = struct;

for iExperiment = 1:numel(results)
    error(iExperiment).theta_zonn = rmse(results(iExperiment).theta_zonn, results(iExperiment).theta);
    error(iExperiment).theta_em_from_zonn = rmse(results(iExperiment).theta_em_from_zonn(end-20:end), results(iExperiment).theta);
    error(iExperiment).theta_em_oracle = rmse(results(iExperiment).theta_em_oracle(end), results(iExperiment).theta);
end

end