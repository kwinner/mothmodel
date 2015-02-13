function plot_mean_observation( mu, sigma, lambda, T_plot, N, alpha, varargin )

DEFAULT_MODEL = 'gaussian'; %gaussian, zonn

parser = inputParser;
addOptional(  parser, 'y',     []);
addOptional(  parser, 'T_obs', []);
addParamValue(parser, 'model', DEFAULT_MODEL);

parser.parse(varargin{:});
y     = parser.Results.y;
T_obs = parser.Results.T_obs;
model = parser.Results.model;

PLOT_TITLE = sprintf('Mean observation, model = %s', model);
XLABEL     = 'Days after May 1st';
YLABEL     = 'Count';

TITLE_SIZE = 14;
LABEL_SIZE = 12;

%T_plot needs to be a superset of T_obs for now
T_plot = union(T_plot, T_obs);

pt_plot = presence_probs(mu, sigma, lambda, T_plot);

figure
hold on
title(PLOT_TITLE, 'FontSize', TITLE_SIZE)
xlabel(XLABEL,    'FontSize', LABEL_SIZE)
ylabel(YLABEL,    'FontSize', LABEL_SIZE)

%compute the mean (same for all models)
mean_plot = alpha * N .* pt_plot;

%under the gaussian model, we can plot the confidence intervals (which need to be drawn first)
if find(ismember({'gaussian', 'GP'}, model)) ~= 0 & ~isempty(T_obs)
	cov_plot = gaussian_cov(pt_plot, lambda, T_plot, N, alpha);

	%group the mean and covariance matrix by obs/unobs
	[T_obs,   i_obs]   = intersect(T_plot, T_obs, 'stable');
	[T_unobs, i_unobs] = setdiff(  T_plot, T_obs, 'stable');

	a = mean_plot(i_obs);
	b = mean_plot(i_unobs);

	A = cov_plot(i_obs,   i_obs);
	B = cov_plot(i_unobs, i_unobs);
	C = cov_plot(i_obs,   i_unobs);

	%compute the conditional distribution for the unobserved times given the obs times
	mean_cond = b' + C' * inv(A) * (y - a)';
	cov_cond  = B  - C' * inv(A) * C;
	var_cond  = cov_cond(eye(size(cov_cond)) == 1);

	%confidence interval width
	conf_cond = sqrt(var_cond) .* 2;

	background = fill([T_unobs'; flipdim(T_unobs', 1)], [mean_cond + conf_cond; flipdim(mean_cond - conf_cond, 1)], [7 7 7]/8);
	background.FaceAlpha = 0.5;
end

%plot the mean
switch model
case {'gaussian', 'GP'}
	plot(T_unobs, mean_cond);
case 'zonn'
	plot(T_plot, mean_plot)
end

%draw the observations, if they exist
if ~isempty(T_obs) && ~isempty(y)
	plot(T_obs, y, '+')
end

end

