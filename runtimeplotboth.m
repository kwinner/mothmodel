%written 10/26/16 by Kevin for NIPS CRC

% close all

% line_width = 1.25;
marker_size = 46;
line_width = 1;
% marker_size = 36;

title_size = 22;
xlabel_size = 20;
ylabel_size = xlabel_size;
legend_size = 18;
axes_size = 16;

x = zeros(nN, nAlpha);
for iN = 1:nN
    for iAlpha = 1:nAlpha
        x(iN, iAlpha) = N_hatVec(iN) * alphaVec(iAlpha);
    end
end

h_fig = figure;
hold on
set(gca, 'FontSize', axes_size)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')


h_pgffa = scatter(x(:), meanRuntimePGFFA(:), marker_size, ...
    'o', ...
    'LineWidth', line_width);

h_fa_or = scatter(x(:), meanRuntimeFAatOracle(:), marker_size, ...
    'd', ...
    'LineWidth', line_width);

h_fa_po = scatter(x(:), meanRuntimeFAatPoiss(:), marker_size, ...
    'x', ...
    'LineWidth', line_width);

title('\Lambda\rho vs Runtime', 'FontSize', title_size)
legend([h_fa_po, h_fa_or, h_pgffa], {'FA - Poiss', 'FA - Oracle', 'PGFFA'}, 'Location', 'northwest', 'FontSize', legend_size)
% legend('show')
xlabel('\Lambda\rho', 'FontSize', xlabel_size)
ylabel('Mean runtime (s)', 'FontSize', ylabel_size)

xtrim = 1.4;

set(h_fig, 'PaperUnits','centimeters');
set(h_fig, 'Units','centimeters');
pos=get(h_fig,'Position');
set(h_fig, 'PaperSize', [(pos(3) - xtrim) pos(4)]);
set(h_fig, 'PaperPositionMode', 'manual');
set(h_fig, 'PaperPosition',[(0 - .25*xtrim) 0 pos(3) pos(4)]);
% set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
print('-dpdf','nalpha_loglog.pdf');