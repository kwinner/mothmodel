function plot_individual_lifetimes( S, Z, varargin )
% PLOT_INDIVIDUAL_LIFETIMES := Draw the lifetimes of every individual stacked above a timeline
% plot_individual_lifetimes( S, Z, ... )
%        (draw all individuals on a timeline)
% plot_individual_lifetimes( S, Z, T, ... )
%        (also draw the times when the population was observed)
%    S                    = vector [N x 1] of individual birth times
%    Z                    = vector [N x 1] of individual death times
%    [T]                  = vector [1 x K] of observation times (sample times)
%    "colors"             = vector [N x 1] of color indices to draw each individual
%                           note: the actual colors are generated from distinguishable_colors.m, these are indices
%    "showAbundance"      = boolean flag whether or not to draw abundance above each observation
%    "observLineWidth"    = line width for vertical time lines
%    "observLineColor"    = color of vertical time lines
%    "abundanceColor"     = color of abundance text
%    "abundanceLabel"     = text to label abundance numbers
%    "abundanceLabelSize" = font size of abundance label
%    "title"              = text to title plot with
%    "titleSize"          = font size of plot title
%    "lifetimeLineWidth"  = line width for individual lines
%    "lifetimeMarkerSize" = marker size of endcaps on lifespan lines
%    "xlabel"             = text for xlabel
%    "xlabelSize"         = font size of xlabel
%    "figureSize"         = size (outerposition) of figure

assert(size(S, 1) == size(Z, 1))

%defaults
DEFAULT_SHOW_ABUNDANCE       = false;
DEFAULT_OBSERV_LINE_WIDTH    = 1;
DEFAULT_OBSERV_LINE_COLOR    = 'r';
DEFAULT_ABUNDANCE_COLOR      = 'r';
DEFAULT_ABUNDANCE_LABEL      = 'Abundance:';
DEFAULT_ABUNDANCE_LABEL_SIZE = 12;

DEFAULT_TITLE                = 'Individual lifetimes';
DEFAULT_TITLE_SIZE           = 14;

DEFAULT_LIFETIME_LINE_WIDTH  = 2;
DEFAULT_LIFETIME_MARKER_SIZE = 2;

DEFAULT_XLABEL               = 'Days after May 1st';
DEFAULT_XLABEL_SIZE          = 12;

DEFAULT_FIGURE_SIZE          = [];

%input parsing
parser = inputParser;
addOptional(  parser, 'T',                  []);
addParamValue(parser, 'colors',             ones(size(S)));
addParamValue(parser, 'showAbundance',      DEFAULT_SHOW_ABUNDANCE);
addParamValue(parser, 'observLineWidth',    DEFAULT_OBSERV_LINE_WIDTH);
addParamValue(parser, 'observLineColor',    DEFAULT_OBSERV_LINE_COLOR);
addParamValue(parser, 'abundanceColor',     DEFAULT_ABUNDANCE_COLOR);
addParamValue(parser, 'abundanceLabel',     DEFAULT_ABUNDANCE_LABEL);
addParamValue(parser, 'abundanceLabelSize', DEFAULT_ABUNDANCE_LABEL_SIZE);
addParamValue(parser, 'title',              DEFAULT_TITLE);
addParamValue(parser, 'titleSize',          DEFAULT_TITLE_SIZE);
addParamValue(parser, 'lifetimeLineWidth',  DEFAULT_LIFETIME_LINE_WIDTH);
addParamValue(parser, 'lifetimeMarkerSize', DEFAULT_LIFETIME_MARKER_SIZE);
addParamValue(parser, 'xlabel',             DEFAULT_XLABEL);
addParamValue(parser, 'xlabelSize',         DEFAULT_XLABEL_SIZE);
addParamValue(parser, 'figureSize',         DEFAULT_FIGURE_SIZE);

parser.parse(varargin{:});
T                  = parser.Results.T;
colors             = parser.Results.colors;
showAbundance      = parser.Results.showAbundance;
observLineWidth    = parser.Results.observLineWidth;
observLineColor    = parser.Results.observLineColor;
abundanceColor     = parser.Results.abundanceColor;
abundanceLabel     = parser.Results.abundanceLabel;
abundanceLabelSize = parser.Results.abundanceLabelSize;
plotTitle          = parser.Results.title;
titleSize          = parser.Results.titleSize;
lifetimeLineWidth  = parser.Results.lifetimeLineWidth;
lifetimeMarkerSize = parser.Results.lifetimeMarkerSize;
xlabelText         = parser.Results.xlabel;
xlabelSize         = parser.Results.xlabelSize;
figureSize         = parser.Results.figureSize;

%readability definitions
N = size(S, 1);
K = size(T, 2);

assert(size(colors, 1) == N)

%stack the individuals (by computing the Y values)
Y = [(1:N); (1:N)];

figure('Size', ); hold on
ax1 = gca;

%draw each group of individuals with the same color
uniqueColors = unique(colors);
plotColors   = distinguishable_colors(numel(uniqueColors));
for icolor = 1:numel(uniqueColors)
	color  = uniqueColors(icolor);
	plotColor = plotColors(icolor, :); %cycle through available colors

	colorGroup = (colors == color); %all individuals with this color

	plot([S(colorGroup), S(colorGroup) + Z(colorGroup)]', Y(:, colorGroup), ...
		 's-', 'LineWidth', lifetimeLineWidth, 'MarkerSize', lifetimeMarkerSize, ...
		 'Color', plotColor);
end

%compute the limits
xlim([min(S) - 1, max(S + Z) + 1])
ylim([0, max(Y(:)) + 1])

%remove the yaxis
set(ax1, 'YTick', [])
set(ax1, 'YColor', get(ax1, 'Color'))

%label stuff
h = xlabel(xlabelText, 'FontSize', xlabelSize);
set(h, 'Units', 'Normalized')
set(h, 'Position', get(h, 'Position') - [0, .025, 0])

titleHandle = title( plotTitle,  'FontSize', titleSize);

%if observation times provided, draw them
if ~isempty(T)
	plot([T;T], repmat(ylim, K, 1)', ...
		 '--', 'LineWidth', observLineWidth, 'Color', observLineColor)

	%replace the xticks with the values in T
	set(gca, 'XTick', T);

	%then, if asked to draw abundance, do so
	if showAbundance
		%compute n from S,Z
		n = abundance(S, Z, T);

		%writing the abundance at the top requires a second axes
		ax2 = axes('Position', get(ax1, 'Position'), ...
			       'XLim', xlim, 'YLim', ylim, ...
			       'Color', 'none', ...
			       'XColor', abundanceColor, ...
			       'YColor', get(ax1, 'Color'), ...
			       'XAxisLocation', 'top');

		%add n as xticks and remove the yticks
		set(ax2, 'XTick', T)
		set(ax2, 'XTickLabel', n)
		set(ax2, 'YTick', [])

		%add a label for abundance
		h = xlabel(abundanceLabel, 'Color', abundanceColor, 'FontSize', abundanceLabelSize);

		%move the label to be in line with the ticks
		set(h, 'Units', 'normalized')
		set(h, 'Position', get(h, 'Position') - [.505, .035, 0])
		set(h, 'HorizontalAlignment', 'right')

		%move the title out of the way
		set(titleHandle, 'Units', 'normalized')
		set(titleHandle, 'Position', get(titleHandle, 'Position') + [0, .04, 0])
	end
end

hold off

end

