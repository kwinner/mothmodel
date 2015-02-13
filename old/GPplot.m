function GPplot (theta, y, params_obs, t_plot)

%throughout this function there are some variables with the _plot suffix and some with _obs
%_obs variables are observed and conditioned on
%_plot variables are being inferred/interpolated for plotting

%clone params over to a new struct for plotting
params_plot   = params_obs;
params_plot.t = t_plot;
T_plot        = numel(params_plot.t);

%compute the pdf for all intervals in the full space
p_plot = ppdf(theta, params_plot);

%combine the pdf entries to get observational likelihoods
pt_plot = arrayfun(@(k) sum(sum(p_plot(1:k, k+1:T_plot+1))), 1:T_plot);

%compute the mean
mean_plot = params_plot.alpha * params_plot.N .* pt_plot;

cov_plot = GPcov(pt_plot, theta, params_plot);

%reorganize the mean and cov according to which variables were observed
[t_obs,i_obs,~] = intersect(params_plot.t, params_obs.t, 'stable');
[t_unobs,i_unobs] = setdiff(params_plot.t, params_obs.t, 'stable');

a = mean_plot(i_obs);
b = mean_plot(i_unobs);

A = cov_plot(i_obs, i_obs);
B = cov_plot(i_unobs, i_unobs);
C = cov_plot(i_obs, i_unobs);

%compute the conditional distribution for the unobserved times given the obs times
mean_cond = b' + C' * inv(A) * (y - a)';
cov_cond  = B - C' * inv(A) * C;
var_cond  = cov_cond(eye(size(cov_cond)) == 1);

%reconstruct plot matrices
conf_cond = sqrt(var_cond) .* 2;

%plot
figure
hold on
fill([t_unobs'; flipdim(t_unobs', 1)], [mean_cond + conf_cond; flipdim(mean_cond - conf_cond, 1)], [7 7 7]/8);
plot(t_unobs, mean_cond);
plot(t_obs, y, '+');

keyboard

end