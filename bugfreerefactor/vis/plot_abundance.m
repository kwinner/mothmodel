function plot_abundance(S, Z, varargin)

DEFAULT_PLOT_TITLE     = 'True abundance over time';
DEFAULT_NEWFIGURE      = true;

ABUNDANCE_LINE_WIDTH  = 1.001;
ABUNDANCE_MARKER_SIZE = 6;
OBS_TIME_LINE_WIDTH   = 1;
OBS_MARKER_SIZE       = 8;

TITLE_SIZE = 14;
LABEL_SIZE = 12;

parser = inputParser;
addOptional(  parser, 'T',              []);
addOptional(  parser, 'y',              []);
addParamValue(parser, 'title',          DEFAULT_PLOT_TITLE);
addParamValue(parser, 'newfigure',      DEFAULT_NEWFIGURE);

parser.parse(varargin{:})
T              = parser.Results.T;
y              = parser.Results.y;
plottitle      = parser.Results.title;
newfigure      = parser.Results.newfigure;

%sometimes we do want to overlay multiple figures
if newfigure; figure; end

N = numel(S);

%put the birth and death times into a single vector
eventtimes = [S; S+Z];
eventsteps = [ones(N,1); -1.*ones(N,1)];

%sort both vectors by eventtimes
[eventtimes, inds] = sort(eventtimes);
eventsteps = eventsteps(inds);

%convert steps into abundance
abundance = cumsum(eventsteps);

stairs([0;eventtimes], [0;abundance], 'LineWidth', ABUNDANCE_LINE_WIDTH);

%plotting bounds
X_LIM = [0,max(S+Z)+1];
Y_LIM = [0,max(abundance) + 1];
xlim(X_LIM);
ylim(Y_LIM);

%label stuff
if newfigure
	xlabel('Days after May 1', 'FontSize', LABEL_SIZE)
	ylabel('Count', 'FontSize', LABEL_SIZE);

	%configure the title
	titlehandle = title(plottitle);
	set(titlehandle, 'FontSize', TITLE_SIZE);
	set(titlehandle, 'Position', get(titlehandle, 'Position') + [0, .3, 0]);
end

%if observation times provided, draw them
if ~isempty(T)
	K = numel(T);

	hold on
	plot([T;T], repmat(Y_LIM,K,1)', '--r', 'LineWidth', OBS_TIME_LINE_WIDTH)

	%replace the ticks with the values in T
	ax = gca;
	set(ax, 'XTick', sort(unique([X_LIM(1):5:X_LIM(2), T])));

	n = arrayfun(@(t) sum(S <= t & S+Z >= t), T);

	%plot markers for these
	plot(T, n, 's', 'MarkerSize', ABUNDANCE_MARKER_SIZE)

	if ~isempty(y)
		plot(T, y, 'og', 'MarkerSize', OBS_MARKER_SIZE, 'MarkerFaceColor', 'g')
	end

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

hold off

end

