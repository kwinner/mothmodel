function plot_true_abundance( S, Z, varargin )
%PLOT_TRUE_ABUNDANCE -- Draw the step function of abundance over time
%params
%required:
%   S         = column vector of birth times
%   Z         = column vector of lifespans
%               S and Z can be 2D matrices or cell arrays of the same dimension
%               in which case each column will be plotted independently
%optional:
%   T         = row vector of observation times
%   y         = row vector of observed counts
%paramvalue:
%   title     = title to draw on the plot
%   lineColor = color of each line
%               if S,Z have k columns, lineColor can be a 2D matrix of size (k,3) or (k,4)
%               note specifying the alpha channel requires HG2 (enabled by default in R2014b)
%   draw_T    = whether to draw lines at the observation times
%   draw_n    = whether to label the true abundance at the observation times
%               if S,Z have multiple columns, only the first values will be labeled
%   draw_y    = whether to put markers at the observed counts

DEFAULT_LINECOLOR     = [];
DEFAULT_DRAW_T        = true;
DEFAULT_DRAW_N        = true;
DEFAULT_DRAW_Y        = false;

PLOT_TITLE            = 'True abundance over time';
XLABEL                = 'Days after May 1st';
YLABEL                = 'Count';

ABUNDANCE_LINE_WIDTH  = 1.001;
ABUNDANCE_MARKER_SIZE = 6;
OBS_TIME_LINE_WIDTH   = 1;
OBS_MARKER_SIZE       = 8;

TITLE_SIZE            = 14;
LABEL_SIZE            = 12;

parser = inputParser;
addOptional(  parser, 'T',              []);
addOptional(  parser, 'y',              []);
addParamValue(parser, 'lineColor',      DEFAULT_LINECOLOR);
addParamValue(parser, 'draw_T',         DEFAULT_DRAW_T);
addParamValue(parser, 'draw_n',         DEFAULT_DRAW_N);
addParamValue(parser, 'draw_y',         DEFAULT_DRAW_Y);

parser.parse(varargin{:})
T              = parser.Results.T;
y              = parser.Results.y;
lineColor      = parser.Results.lineColor;
draw_T         = parser.Results.draw_T;
draw_n         = parser.Results.draw_n;
draw_y         = parser.Results.draw_y;

%if S, Z are not cell arrays, wrap them in cell arrays
%using a 2D matrix would be cleaner, but would require every S,Z have the same number of individuals
S = num2cell(S, 1);
Z = num2cell(Z, 1);

nSamples = numel(S);

N = cellfun(@rows, S);

%put the birth and death times into a single vector
eventtimes = cellfun(@(s,z) [s;s+z], S, Z, 'UniformOutput', false);
eventsteps = arrayfun(@(N_i) [ones(N_i,1); -1.*ones(N_i,1)], N, 'UniformOutput', false);

%sort both vectors by eventtimes
[eventtimes, inds] = cellfun(@sort,eventtimes, 'UniformOutput', false);
eventsteps = cellfun(@(eventsteps_col, inds_col) eventsteps_col(inds_col), eventsteps, inds, 'UniformOutput', false);

%convert steps into abundance
abundance = cellfun(@cumsum, eventsteps, 'UniformOutput', false);

figure
hold on
for iSample = 1:nSamples
	[stairs_x,stairs_y] = stairs([0;eventtimes{iSample}], [0;abundance{iSample}]);
	nLine = plot(stairs_x, stairs_y, 'LineWidth', ABUNDANCE_LINE_WIDTH);

	%with HG2 (standard in R2014b, activatible w/ cmd line param -hgVersion 2 in earlier versions) we can draw transparent lines
	%note: this is an undocumented feature and may break in the future
	if ~isempty(lineColor)
		nLine.Color = lineColor(mod(iSample-1, rows(lineColor))+1,:);
	end
end

%plotting bounds
X_LIM = [0,max(cellfun(@(s,z) max(s+z), S, Z))+1];
Y_LIM = [0,max(cellfun(@max, abundance)) + 1];
xlim(X_LIM);
ylim(Y_LIM);

%label stuff
xlabel(XLABEL, 'FontSize', LABEL_SIZE)
ylabel(YLABEL, 'FontSize', LABEL_SIZE);

%configure the title
titlehandle = title(PLOT_TITLE);
set(titlehandle, 'FontSize', TITLE_SIZE);
set(titlehandle, 'Position', get(titlehandle, 'Position') + [0, .3, 0]);

%if observation times provided, draw them
if ~isempty(T)
	K = numel(T);

	if draw_T
		plot([T;T], repmat(Y_LIM,K,1)', '--r', 'LineWidth', OBS_TIME_LINE_WIDTH)

		%replace the ticks with the values in T
		ax = gca;
		set(ax, 'XTick', sort(unique([X_LIM(1):5:X_LIM(2), T])));
	end

	if draw_y && ~isempty(y)
		plot(T, y, 'og', 'MarkerSize', OBS_MARKER_SIZE, 'MarkerFaceColor', 'g')
	end

	if draw_n
		n = arrayfun(@(t) sum(S{1} <= t & S{1}+Z{1} >= t), T);

		%plot markers for these
		plot(T, n, 's', 'MarkerSize', ABUNDANCE_MARKER_SIZE)

		%clone the existing axes
		ax2 = axes('Position', get(ax, 'Position'), 'Color', 'none');
		set(ax2, 'XAxisLocation', 'top');
		set(ax2, 'XLim', get(ax, 'XLim'));

		%draw n as labels over T
		set(ax2, 'XTick', T);
		set(ax2, 'XTickLabel', n);
		set(ax2, 'YTick', []);

		set(ax2, 'XColor', 'r');

		%label it?
		xlabel(ax2,'Abundance:','color','r')
		xlabelhandle = get(ax2, 'XLabel');
		set(xlabelhandle, 'Position', [0, 1.01, 0])
		set(xlabelhandle, 'horizontalAlignment', 'left')
	end
end

end

