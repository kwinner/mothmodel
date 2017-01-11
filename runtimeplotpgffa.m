%written 10/26/16 by Kevin for NIPS CRC

% close all

% line_width = 1.25;
% marker_size = 60;
line_width = 1;
marker_size = 36;

title_size = 22;
xlabel_size = 20;
ylabel_size = xlabel_size;
legend_size = 20;
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

scatter(x(:), meanRuntimePGFFA(:), marker_size, ...
    'LineWidth', line_width)

title('\Lambda\rho vs Runtime of PGFFA', 'FontSize', title_size)
legend({'PGFFA'}, 'Location', 'northwest', 'FontSize', legend_size)
xlabel('\Lambda\rho', 'FontSize', xlabel_size)
ylabel('Mean runtime (s)', 'FontSize', ylabel_size)

xtrim = 3;

set(h_fig, 'PaperUnits','centimeters');
set(h_fig, 'Units','centimeters');
pos=get(h_fig,'Position');
set(h_fig, 'PaperSize', [(pos(3) - xtrim) pos(4)]);
set(h_fig, 'PaperPositionMode', 'manual');
set(h_fig, 'PaperPosition',[(0 - .25*xtrim) 0 (pos(3)-xtrim/2) pos(4)]);
% set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
print('-dpdf','nalpha_pgffa.pdf');
