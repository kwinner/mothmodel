function plot_true_abundance( S, Z, varargin )
% PLOT_TRUE_ABUNDANCE := Draw the step function of abundance over time
% plot_true_abundance( S, Z, ... )
% plot_true_abundance( S, Z, T, ... )
% plot_true_abundance( S, Z, T, y, ... )
%
% INPUTS
% required:
%    S              = vector [N x 1] of individual birth times
%    Z              = vector [N x 1] of individual lifespans
%                     S and Z can be 2D matrices of size [N x nSamples]
%                     in which case each column will be plotted independently
%                     S and Z can also be cell arrays [1 x nSamples]
%                     in which case each vector will be plotted independently (and can have different N)
% optional:
%    [T]            = vector [1 x K] of observation times
%    [y]            = vector [1 x K] of observed counts
% paramvalue:
%    'lineColor'    = color of each line
%                     if S,Z have k columns, lineColor can be a 2D matrix of size (k,3) or (k,4)
%                     note: specifying the alpha channel requires HG2 (enabled by default in R2014b)
%    'lineOpacity'  = opacity of each line (doesn't work in all versions)
%                     note: if the colors of lineColor specify the alpha channel, this is overrided
%    'lineWidth'    = line width for abundance line
%    'draw_T'       = whether to draw lines at the observation times
%    'TWidth'       = line width of observation time lines
%    'TColor'       = color of observation time lines
%    'TStyle'       = style of observation time lines
%    'draw_n'       = whether to label the true abundance at the observation times
%                     if S,Z have multiple columns, only the first values will be labeled
%    'nMarkerStyle' = marker style of abundance
%    'nMarkerColor' = marker color of abundance
%    'nMarkerSize'  = marker size of abundance
%    'draw_y'       = whether to put markers at the observed counts
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

%defaults
DEFAULT_LINECOLOR     = [];
DEFAULT_LINEOPACITY   = [];
DEFAULT_LINEWIDTH     = 1.001;

DEFAULT_DRAW_T        = true;
DEFAULT_TWIDTH        = 1;
DEFAULT_TCOLOR        = 'r';
DEFAULT_TSTYLE        = '--';

DEFAULT_DRAW_N        = true;
DEFAULT_NMARKERSTYLE  = 's';
DEFAULT_NMARKERCOLOR  = 'r';
DEFAULT_NMARKERSIZE   = 6;

%if y is provided, it will always be drawn (no reason to include it otherwise)
DEFAULT_YMARKERSTYLE  = '+';
DEFAULT_YMARKERCOLOR  = 'b';
DEFAULT_YMARKERSIZE   = 6;

DEFAULT_TITLE         = 'True abundance over time';
DEFAULT_TITLESIZE     = 14;

DEFAULT_XLABEL        = 'Days after May 1st';
DEFAULT_XLABELSIZE    = 12;

DEFAULT_YLABEL        = 'Abundance';
DEFAULT_YLABELSIZE    = 12;

DEFAULT_TICKS         = [];
DEFAULT_TICK_PATTERN  = 'add';
POSSIBLE_TICK_PATTERN = {'add','replace','none'};

DEFAULT_FIGURE_SIZE   = [];

parser = inputParser;
addOptional(  parser, 'T',              []);
addOptional(  parser, 'y',              []);
addParamValue(parser, 'lineColor',    DEFAULT_LINECOLOR);
addParamValue(parser, 'lineOpacity',  DEFAULT_LINEOPACITY);
addParamValue(parser, 'lineWidth',    DEFAULT_LINEWIDTH);
addParamValue(parser, 'draw_T',       DEFAULT_DRAW_T);
addParamValue(parser, 'TWidth',       DEFAULT_TWIDTH);
addParamValue(parser, 'TColor',       DEFAULT_TCOLOR);
addParamValue(parser, 'TStyle',       DEFAULT_TSTYLE);
addParamValue(parser, 'draw_n',       DEFAULT_DRAW_N);
addParamValue(parser, 'nMarkerStyle', DEFAULT_NMARKERSTYLE);
addParamValue(parser, 'nMarkerColor', DEFAULT_NMARKERCOLOR);
addParamValue(parser, 'nMarkerSize',  DEFAULT_NMARKERSIZE);
addParamValue(parser, 'yMarkerStyle', DEFAULT_YMARKERSTYLE);
addParamValue(parser, 'yMarkerColor', DEFAULT_YMARKERCOLOR);
addParamValue(parser, 'yMarkerSize',  DEFAULT_YMARKERSIZE);
addParamValue(parser, 'title',        DEFAULT_TITLE);
addParamValue(parser, 'titleSize',    DEFAULT_TITLESIZE);
addParamValue(parser, 'xlabel',       DEFAULT_XLABEL);
addParamValue(parser, 'xlabelSize',   DEFAULT_XLABELSIZE);
addParamValue(parser, 'ylabel',       DEFAULT_YLABEL);
addParamValue(parser, 'ylabelSize',   DEFAULT_YLABELSIZE);
addParamValue(parser, 'ticks',        DEFAULT_TICKS);
addParamValue(parser, 'tickPattern',  DEFAULT_TICK_PATTERN);
addParamValue(parser, 'figureSize',   DEFAULT_FIGURE_SIZE);

parser.parse(varargin{:})
T            = parser.Results.T;
y            = parser.Results.y;
lineColor    = parser.Results.lineColor;
lineOpacity  = parser.Results.lineOpacity;
lineWidth    = parser.Results.lineWidth;
draw_T       = parser.Results.draw_T;
TWidth       = parser.Results.TWidth;
TColor       = parser.Results.TColor;
TStyle       = parser.Results.TStyle;
draw_n       = parser.Results.draw_n;
nMarkerStyle = parser.Results.nMarkerStyle;
nMarkerColor = parser.Results.nMarkerColor;
nMarkerSize  = parser.Results.nMarkerSize;
yMarkerStyle = parser.Results.yMarkerStyle;
yMarkerColor = parser.Results.yMarkerColor;
yMarkerSize  = parser.Results.yMarkerSize;
titleText    = parser.Results.title;
titleSize    = parser.Results.titleSize;
xlabelText   = parser.Results.xlabel;
xlabelSize   = parser.Results.xlabelSize;
ylabelText   = parser.Results.ylabel;
ylabelSize   = parser.Results.ylabelSize;
ticks        = parser.Results.ticks;
tickPattern  = parser.Results.tickPattern;
figureSize   = parser.Results.figureSize;

%if S, Z are not cell arrays, wrap them in cell arrays
%using a 2D matrix would be cleaner, but would require every S,Z have the same number of individuals
assert(iscell(S) == iscell(Z));
if ~iscell(S)
	S = num2cell(S, 1);
	Z = num2cell(Z, 1);
end

%N is a vector of the number of individuals in each sample
%nSamples counts how many unique samples there are (each will be its own plot)
assert(all(cellfun(@numel, S) == cellfun(@numel, Z)));
N = cellfun(@numel, S);
nSamples = numel(N);

%put the birth and death times into a single vector [2N x 1] for each sample
eventtimes = cellfun(@(s,z) ...
	                 [s;s+z] ...
	                 , S, Z, 'UniformOutput', false);
eventsteps = arrayfun(@(N_i) ...
                      [ones(N_i,1); -1.*ones(N_i,1)] ...
                      , N, 'UniformOutput', false); %is the event a birth or death?

%sort both vectors by eventtimes
[eventtimes, inds] = cellfun(@sort, eventtimes, 'UniformOutput', false);
eventsteps = cellfun(@(eventsteps_sample, inds_sample) ...
	                 eventsteps_sample(inds_sample) ...
	                 , eventsteps, inds, 'UniformOutput', false);

%convert steps into abundance
nTrue = cellfun(@cumsum, eventsteps, 'UniformOutput', false);

%append an initial t0 event with 0 abundance to each sample
eventtimes = cellfun(@(eventtimes_i) ...
	                 [eventtimes_i(1); eventtimes_i] ...
	                 , eventtimes, 'UniformOutput', false);
nTrue      = cellfun(@(nTrue_i) ...
	                 [0; nTrue_i] ...
	                 , nTrue, 'UniformOutput', false);

figure
hold on
for iSample = 1:nSamples
	[stairs_x,stairs_y] = stairs(eventtimes{iSample}, nTrue{iSample});
	nLine = plot(stairs_x, stairs_y, ...
		         'LineWidth', lineWidth);

	%with HG2 (standard in R2014b, activatible w/ cmd line param -hgVersion 2 in earlier versions) we can draw transparent lines
	%note: this is an undocumented feature and may break in the future
	if ~isempty(lineColor)
		%note: line colors cycle through those provided in lineColor, if necessary
		nLine.Color = lineColor(mod(iSample-1, size(lineColor,1))+1,:);
	elseif (isempty(lineColor) || size(lineColor, 2)) ~= 4 && ~isempty(lineOpacity)
		%if opacity is unset by lineColor, and opacity is provided...
		nLine.Color(4) = lineOpacity;
	end
end

%plotting bounds
X_LIM = [min(cellfun(@min, S)) - 1, ...
         max(cellfun(@(s,z) max(s+z), S, Z)) + 1];
Y_LIM = [0,max(cellfun(@max, nTrue)) + 1];
xlim(X_LIM);
ylim(Y_LIM);

%label stuff
title( titleText,  'FontSize', titleSize)
xlabel(xlabelText, 'FontSize', xlabelSize)
ylabel(ylabelText, 'FontSize', ylabelSize)

%overwrite the xticks
if ~isempty(ticks)
	set(gca, 'XTick', ticks);
end

%if observation times provided, we can draw T, n or y
if ~isempty(T)
	K = numel(T);

	%if y is provided, it will always be drawn
	draw_y = ~isempty(y);

	%update xticks
	if any([draw_T, draw_n, draw_y])
		switch validatestring(tickPattern, POSSIBLE_TICK_PATTERN)
		case 'none'
		case 'add'
			%add T to the ticks with the values in T
			set(gca, 'XTick', sort(unique([get(gca, 'XTick'), T])));
		case 'replace'
			%replace ticks with the values in T
			set(gca, 'XTick', sort(T));
		end
	end

	%draw lines at observation times
	if draw_T
		%draw vertical lines at T
		plot([T;T], repmat(Y_LIM,K,1)', ...
			 TStyle, 'LineWidth', TWidth, 'Color', TColor)
	end

	%draw markers at the abundance points
	if draw_n
		n = abundance(S{1}, Z{1}, T);

		plot(T, n, ...
			 nMarkerStyle, 'MarkerSize', nMarkerSize, 'MarkerEdgeColor', nMarkerColor)
	end

	%draw markers at the observed values
	if draw_y
		plot(T, y, ...
			 yMarkerStyle, 'MarkerSize', yMarkerSize, 'MarkerEdgeColor', yMarkerColor)
	end
end

end

