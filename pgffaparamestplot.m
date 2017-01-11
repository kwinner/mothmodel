%written 10/26/16 by Kevin for NIPS CRC

close all

line_width = 2.5;
marker_size = 10;

title_size = 22;
xlabel_size = 20;
ylabel_size = xlabel_size;
legend_size = 18;
axes_size = 14;

fill = false;
trunc_error = false;

x = 10:10:150;

h_fig = figure
hold on
set(gca, 'FontSize', axes_size)

h_pgffa = errorbar(x, mean(N_hatPGFFA_all, 2), N_hatPGFFA_error, ...
    'o-', ...
    'LineWidth', line_width, ...
    'MarkerSize', marker_size);
if fill
    set(h_pgffa, 'MarkerFaceColor', get(h_pgffa, 'Color'));
end

h_trunc = 0;
if trunc_error
    h_trunc = errorbar(x, mean(N_hatTrunc_all, 2), N_hatTrunc_error, ...
    's-', ...
    'LineWidth', line_width, ...
    'MarkerSize', marker_size);
else
    h_trunc = plot(x, mean(N_hatTrunc_all, 2), ...
    'd-', ...
    'LineWidth', line_width, ...
    'MarkerSize', marker_size);
end

if fill
    set(h_trunc, 'MarkerFaceColor', get(h_trunc, 'Color'));
end

h_true = plot(x, x, ...
    '--', ...
    'LineWidth', line_width);

xlim([5,155])
ylim([0,205])
set(gca, 'XTick', 10:20:150)

title('\lambda Recovery', 'FontSize', title_size)
legend([h_trunc, h_pgffa, h_true], {'Trunc', 'PGFFA', 'True \lambda'}, 'Location', 'northwest', 'FontSize', legend_size)
xlabel('\lambda', 'FontSize', xlabel_size)
ylabel('$\hat{\lambda}$','Interpreter', 'latex', 'FontSize', ylabel_size)

xtrim = 1.4;

set(h_fig, 'PaperUnits','centimeters');
set(h_fig, 'Units','centimeters');
pos=get(h_fig,'Position');
set(h_fig, 'PaperSize', [(pos(3) - xtrim) pos(4)]);
set(h_fig, 'PaperPositionMode', 'manual');
set(h_fig, 'PaperPosition',[(0 - .25*xtrim) 0 pos(3) pos(4)]);
% set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
print('-dpdf','paramest.pdf');
