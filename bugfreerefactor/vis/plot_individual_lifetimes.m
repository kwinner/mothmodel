function plot_individual_lifetimes(S, Z, varargin)

DEFAULT_SHOW_ABUNDANCE = true;
DEFAULT_PLOT_TITLE     = 'Individual lifetimes';

PLOT_COLORS            = 'brgcmy';

INDIVIDUAL_LINE_WIDTH  = 2;
INDIVIDUAL_MARKER_SIZE = 2;
OBS_TIME_LINE_WIDTH    = 1;

TITLE_SIZE             = 14;
XLABEL_SIZE            = 12;

N = numel(S);

parser = inputParser;
addOptional(  parser, 'T',              []);
addParamValue(parser, 'show_abundance', DEFAULT_SHOW_ABUNDANCE);
addParamValue(parser, 'title',          DEFAULT_PLOT_TITLE);
addParamValue(parser, 'colors',         ones(N,1));

parser.parse(varargin{:});
T              = parser.Results.T;
show_abundance = parser.Results.show_abundance;
plottitle      = parser.Results.title;
colors         = parser.Results.colors;

%give each individual a number
I = [(1:N);(1:N)];

%draw a bar for each individual
uniquecolors = unique(colors);
figure; hold on
for icolor = 1:numel(uniquecolors)
	color = uniquecolors(icolor);
	%if there was only one color specified, use it to draw everything
	if numel(uniquecolors) == 1
		plot([S, S+Z]', I, [PLOT_COLORS(color), 's-'], 'LineWidth', INDIVIDUAL_LINE_WIDTH, 'MarkerSize', INDIVIDUAL_MARKER_SIZE)
	else
		colorsubset = colors == color;
		plot([S(colorsubset),S(colorsubset)+Z(colorsubset)]', I(:,colorsubset), [PLOT_COLORS(color), 's-'], 'LineWidth', INDIVIDUAL_LINE_WIDTH, 'MarkerSize', INDIVIDUAL_MARKER_SIZE)
	end
end

%plotting bounds
X_LIM = [0,max(S+Z)+1];
Y_LIM = [0,N+1];
% xlim(X_LIM)
% ylim(Y_LIM)
xlim([0 25.8039])
ylim([0,13])

%remove the yticks
ax = gca;
set(ax, 'YTick', []);

%label stuff
xlabel('Days after May 1', 'FontSize', XLABEL_SIZE)

%configure the title
titlehandle = title(plottitle);
set(titlehandle, 'FontSize', TITLE_SIZE);
set(titlehandle, 'Position', get(titlehandle, 'Position') + [0, .5, 0]);

%if observation times provided, draw them
if ~isempty(T)
	K = numel(T);

	plot([T;T], repmat(Y_LIM,K,1)', '--r', 'LineWidth', OBS_TIME_LINE_WIDTH)

	%replace the ticks with the values in T
	set(ax, 'XTick', sort(unique([X_LIM(1):5:X_LIM(2), T])));

	%if abundancies are provided, put them on top of the plot
	if show_abundance
		n = arrayfun(@(t) sum(S <= t & S+Z >= t), T);

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

hold off

end

