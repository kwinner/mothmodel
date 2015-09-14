function plot_mean_observations( arrivalDistn, serviceDistn, N, alpha, T, varargin )
% PLOT_MEAN_OBSERVATIONS := Draw the expected abundance over time
% plot_mean_observations( arrivalDistn, serviceDistn, N_hat, alpha, T, ... )
% plot_mean_observations( arrivalDistn, serviceDistn, N_hat, alpha, T, y, T_obs, ... )
%
% INPUTS
% required:
%    arrivalDistn   = a distribution object for the birth process (typically normal)
%    serviceDistn   = a distribution object for the death process (typically exponential)
%    N_hat          = mean (or actual) superpopulation size (positive int)
%    alpha          = detection probability (probability)
%    T              = vector [1 x K] of observation times (sample times)
% optional:
%    [y]            = vector [1 x K] of observed counts
%    [T_obs]        = vector [1 x K] of observation times
% paramvalue:
%    'lineStyle'    = abundance line style
%    'lineColor'    = abundance line color
%    'lineWidth'    = abundance line width
%    'yMarkerStyle' = marker style of observations
%    'yMarkerColor' = marker color of observations
%    'yMarkerSize'  = marker size of observations
%    'title'        = title to draw on the plot
%    'titleSize'    = font size of plot title
%    'xlabel'       = text for xlabel
%    'xlabelSize'   = font size of xlabel
%    'ylabel'       = text for ylabel
%    'ylabelSize'   = font size of ylabel
%    'figureSize'   = size (outerposition) of figure

DEFAULT_LINESTYLE     = '-';
DEFAULT_LINECOLOR     = [];
DEFAULT_LINEWIDTH     = 4;

DEFAULT_YMARKERSTYLE  = '+';
DEFAULT_YMARKERCOLOR  = 'b';
DEFAULT_YMARKERSIZE   = 6;

DEFAULT_TITLE         = 'Expected abundance over time';
DEFAULT_TITLESIZE     = 14;

DEFAULT_XLABEL        = 'Days after May 1st';
DEFAULT_XLABELSIZE    = 12;

DEFAULT_YLABEL        = 'Abundance';
DEFAULT_YLABELSIZE    = 12;

DEFAULT_FIGURE_SIZE   = [];

parser = inputParser;
addOptional(  parser, 'y',            []);
addOptional(  parser, 'T_obs',        []);
addParamValue(parser, 'lineStyle',    DEFAULT_LINESTYLE);
addParamValue(parser, 'lineColor',    DEFAULT_LINECOLOR);
addParamValue(parser, 'lineWidth',    DEFAULT_LINEWIDTH);
addParamValue(parser, 'yMarkerStyle', DEFAULT_YMARKERSTYLE);
addParamValue(parser, 'yMarkerColor', DEFAULT_YMARKERCOLOR);
addParamValue(parser, 'yMarkerSize',  DEFAULT_YMARKERSIZE);
addParamValue(parser, 'title',        DEFAULT_TITLE);
addParamValue(parser, 'titleSize',    DEFAULT_TITLESIZE);
addParamValue(parser, 'xlabel',       DEFAULT_XLABEL);
addParamValue(parser, 'xlabelSize',   DEFAULT_XLABELSIZE);
addParamValue(parser, 'ylabel',       DEFAULT_YLABEL);
addParamValue(parser, 'ylabelSize',   DEFAULT_YLABELSIZE);
addParamValue(parser, 'figureSize',   DEFAULT_FIGURE_SIZE);

parser.parse(varargin{:});
y            = parser.Results.y;
T_obs        = parser.Results.T_obs;
lineStyle    = parser.Results.lineStyle;
lineColor    = parser.Results.lineColor;
lineWidth    = parser.Results.lineWidth;
yMarkerStyle = parser.Results.yMarkerStyle;
yMarkerColor = parser.Results.yMarkerColor;
yMarkerSize  = parser.Results.yMarkerSize;
titleText    = parser.Results.title;
titleSize    = parser.Results.titleSize;
xlabelText   = parser.Results.xlabel;
xlabelSize   = parser.Results.xlabelSize;
ylabelText   = parser.Results.ylabel;
ylabelSize   = parser.Results.ylabelSize;
figureSize   = parser.Results.figureSize;

%compute the mean abundance at all times
pt = presence_prob(arrivalDistn, serviceDistn, T);
mean_n = alpha * N .* pt;

figure
hold on
line = plot(T, mean_n, ...
            lineStyle, 'LineWidth', lineWidth);
if ~isempty(lineColor)
	line.Color = lineColor;
end

title( titleText,  'FontSize', titleSize)
xlabel(xlabelText, 'FontSize', xlabelSize)
ylabel(ylabelText, 'FontSize', ylabelSize)

%draw the observations, if they exist
if ~isempty(T_obs) && ~isempty(y)
	plot(T_obs, y, ...
		 yMarkerStyle, 'MarkerSize', yMarkerSize, 'MarkerEdgeColor', yMarkerColor)
end

%under the gaussian model, we can plot the confidence intervals (which need to be drawn first)
% if find(ismember({'gaussian', 'GP'}, model)) ~= 0 & ~isempty(T_obs)
% 	cov_plot = gaussian_cov(pt_plot, lambda, T_plot, N, alpha);

% 	%group the mean and covariance matrix by obs/unobs
% 	[T_obs,   i_obs]   = intersect(T_plot, T_obs, 'stable');
% 	[T_unobs, i_unobs] = setdiff(  T_plot, T_obs, 'stable');

% 	a = mean_plot(i_obs);
% 	b = mean_plot(i_unobs);

% 	A = cov_plot(i_obs,   i_obs);
% 	B = cov_plot(i_unobs, i_unobs);
% 	C = cov_plot(i_obs,   i_unobs);

% 	%compute the conditional distribution for the unobserved times given the obs times
% 	mean_cond = b' + C' * inv(A) * (y - a)';
% 	cov_cond  = B  - C' * inv(A) * C;
% 	var_cond  = cov_cond(eye(size(cov_cond)) == 1);

% 	%confidence interval width
% 	conf_cond = sqrt(var_cond) .* 2;

% 	background = fill([T_unobs'; flipdim(T_unobs', 1)], [mean_cond + conf_cond; flipdim(mean_cond - conf_cond, 1)], [7 7 7]/8);
% 	background.FaceAlpha = 0.5;
% end

end

