%written 1/11/16
function icml_runtime_pgffa_vs_trunc_plot(N_hatVec, ...
                                          alphaVec, ...
                                          runtimeUTPPGFFA, ...
                                          runtimePGFFA, ...
                                          runtimeFAatOracle, ...
                                          runtimeFAatPoiss, ...
                                          runtimeFAatNegBin, ...
                                          resultDir)
    nN = numel(N_hatVec);
    nAlpha = numel(alphaVec);

    LINE_LINE_WIDTH  = 4;
    LINE_LABEL_FONT  = 20;
    LINE_TITLE_FONT  = 22;
    LINE_LEGEND_FONT = 12;
    LINE_FIG_SIZE    = [0 0 8 4];
    
    SCATTER_MARKER_SIZE  = 46;
    SCATTER_LINE_WIDTH   = 1;
    SCATTER_TITLE_FONT   = LINE_TITLE_FONT;
    SCATTER_LABEL_FONT   = LINE_LABEL_FONT;
    SCATTER_LEGEND_FONT  = 18;
    SCATTER_AXES_SIZE    = 16;
    SCATTER_LOGLOG_XTRIM = 1.4;
    SCATTER_LINEAR_XTRIM = 3.0;
    
    %take the average over all iterations (for plotting)
    meanRuntimeUTPPGFFA   = mean(runtimeUTPPGFFA,   3);
    meanRuntimePGFFA      = mean(runtimePGFFA,      3);
    meanRuntimeFAatOracle = mean(runtimeFAatOracle, 3);
    meanRuntimeFAatPoiss  = mean(runtimeFAatPoiss,  3);
    meanRuntimeFAatNegBin = mean(runtimeFAatNegBin, 3);
    
    %%%LINE PLOTS%%%
    
    %plot the alpha vs runtime plots (one for each val of N_hat)
    for iN = 1:nN
        N_hat = N_hatVec(iN);

        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(alphaVec, meanRuntimeFAatPoiss(iN, :),  'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(5,:))
        plot(alphaVec, meanRuntimeFAatNegBin(iN, :), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(4,:))
        plot(alphaVec, meanRuntimeFAatOracle(iN, :), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(3,:))
        plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(2,:))
        plot(alphaVec, meanRuntimeUTPPGFFA(iN, :),   'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(2,:))

        % title(sprintf('$$\\alpha$$ vs runtime of FA and PGFFA, $$\\hat{N} = %d$$', N_hat), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('alpha vs runtime of all methods, N = %d', N_hat), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('alpha',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'FA - Poiss'; ...
                'FA - NegBin'; ...
                'FA - Oracle'; ...
                'PGFFA'; ...
                'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('AllMethods_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
        close(h_fig)

        %repeat for oracle vs PGFFA only
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(alphaVec, meanRuntimeFAatOracle(iN, :), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(3,:))
        plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(2,:))
        plot(alphaVec, meanRuntimeUTPPGFFA(iN, :),   'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\alpha$$ vs runtime of FA and PGFFA, $$\\hat{N} = %d$$', N_hat), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('alpha vs runtime of all methods, N = %d', N_hat), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('alpha',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'FA - Oracle'; ...
                'PGFFA'; ...
                'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('Oracle_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
        close(h_fig)

        %repeat for PGFFA only
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(alphaVec, meanRuntimePGFFA(iN, :),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\alpha$$ vs runtime of PGFFA, $$\\hat{N} = %d$$', N_hat), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('alpha vs runtime of PGFFA, N = %d', N_hat), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('alpha',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'PGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('PGFFA_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
        close(h_fig)

        %repeat for UTPPGFFA only
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(alphaVec, meanRuntimeUTPPGFFA(iN, :),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\alpha$$ vs runtime of PGFFA, $$\\hat{N} = %d$$', N_hat), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\alpha$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('alpha vs runtime of UTPPGFFA, N = %d', N_hat), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('alpha',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('UTPPGFFA_AlphavsRuntime_N%d', N_hat)), ...
              '-depsc')
        close(h_fig)
    end

    %plot the N_hat vs runtime plots (one for each val of alpha)
    for iAlpha = 1:nAlpha
        alpha = alphaVec(iAlpha);

        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(N_hatVec, meanRuntimeFAatPoiss(:, iAlpha),  'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(5,:))
        plot(N_hatVec, meanRuntimeFAatNegBin(:, iAlpha), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(4,:))
        plot(N_hatVec, meanRuntimeFAatOracle(:, iAlpha), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(3,:))
        plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(2,:))
        plot(N_hatVec, meanRuntimeUTPPGFFA(:, iAlpha),   'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\hat{N}$$ vs runtime of FA and PGFFA, $$\\alpha = %0.2f$$', alpha), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('N vs runtime of all methods, alpha = %d', alpha), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('N',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'FA - Poiss'; ...
                'FA - NegBin'; ...
                'FA - Oracle'; ...
                'PGFFA'; ...
                'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('AllMethods_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
        close(h_fig)

        %repeat for oracle vs pgffa vs UTPPGFFA
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(N_hatVec, meanRuntimeFAatOracle(:, iAlpha), 'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(3,:))
        plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(2,:))
        plot(N_hatVec, meanRuntimeUTPPGFFA(:, iAlpha),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\hat{N}$$ vs runtime of FA and PGFFA, $$\\alpha = %0.2f$$', alpha), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)
        title(sprintf('N vs runtime of all methods, alpha = %d', alpha), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('N',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'FA - Oracle'; ...
                'PGFFA'; ...
                'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('Oracle_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
        close(h_fig)

        %repeat for pgffa only
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(N_hatVec, meanRuntimePGFFA(:, iAlpha),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\hat{N}$$ vs runtime of PGFFA, $$\\alpha = %0.2f$$', alpha), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)

        title(sprintf('N vs runtime of PGFFA, alpha = %d', alpha), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('N',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'PGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('PGFFA_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
        close(h_fig)
        
        %repeat for utppgffa only
        % hFig = figure('Visible', 'off');
        % hFig.PaperUnits = 'inches';
        % hFig.PaperPosition = FIGSIZE;
        h_fig = figure('Visible', 'off', 'Units', 'inches', 'Position', LINE_FIG_SIZE);

        hold on
        colOrd = get(gca, 'ColorOrder');

        plot(N_hatVec, meanRuntimeUTPPGFFA(:, iAlpha),      'LineWidth', LINE_LINE_WIDTH, 'Color', colOrd(1,:))

        % title(sprintf('$$\\hat{N}$$ vs runtime of PGFFA, $$\\alpha = %0.2f$$', alpha), ...
        %       'FontSize', TITLEFONT, 'Interpreter', 'Latex')
        % xlabel('$$\hat{N}$$',       'FontSize', LABELFONT, 'Interpreter', 'Latex')
        % ylabel('mean runtime (s)', 'FontSize', LABELFONT)

        title(sprintf('N vs runtime of UTPPGFFA, alpha = %d', alpha), ...
              'FontSize', LINE_TITLE_FONT)
        xlabel('N',       'FontSize', LINE_LABEL_FONT)
        ylabel('mean runtime (s)', 'FontSize', LINE_LABEL_FONT)

        legend({'UTPPGFFA'}, ...
                'FontSize',LINE_LEGEND_FONT,'Location','Northwest')
        hold off

        print(h_fig, ...
              fullfile(resultDir, ...
                       sprintf('UTPPGFFA_NhatvsRuntime_alpha%0.2f.eps', alpha)), ...
              '-depsc')
        close(h_fig)
    end
    
    %%%SCATTER PLOTS%%%
    %compute new xaxis values from N*alpha
    x = zeros(nN, nAlpha);
    for iN = 1:nN
        for iAlpha = 1:nAlpha
            x(iN, iAlpha) = N_hatVec(iN) * alphaVec(iAlpha);
        end
    end
    
    %scatter all methods
    h_fig = figure('Visible', 'off');
    hold on
    set(gca, 'FontSize', SCATTER_AXES_SIZE)
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    
    h_utppgffa = scatter(x(:), meanRuntimeUTPPGFFA(:), SCATTER_MARKER_SIZE, ...
                         'o', ...
                         'LineWidth', SCATTER_LINE_WIDTH);
    
    h_pgffa = scatter(x(:), meanRuntimePGFFA(:), SCATTER_MARKER_SIZE, ...
                      'o', ...
                      'LineWidth', SCATTER_LINE_WIDTH);

    h_fa_or = scatter(x(:), meanRuntimeFAatOracle(:), SCATTER_MARKER_SIZE, ...
                      'd', ...
                      'LineWidth', SCATTER_LINE_WIDTH);

    h_fa_po = scatter(x(:), meanRuntimeFAatPoiss(:), SCATTER_MARKER_SIZE, ...
                      'x', ...
                      'LineWidth', SCATTER_LINE_WIDTH);

    title('\Lambda\rho vs Runtime', 'FontSize', SCATTER_TITLE_FONT)
    legend([h_fa_po, h_fa_or, h_pgffa, h_utppgffa], {'FA - Poiss', 'FA - Oracle', 'PGFFA', 'UTPPGFFA'}, 'Location', 'northwest', 'FontSize', SCATTER_LEGEND_FONT)
    % legend('show')
    xlabel('\Lambda\rho', 'FontSize', SCATTER_LABEL_FONT)
    ylabel('Mean runtime (s)', 'FontSize', SCATTER_LABEL_FONT)

    set(h_fig, 'PaperUnits','centimeters');
    set(h_fig, 'Units','centimeters');
    pos=get(h_fig,'Position');
    set(h_fig, 'PaperSize', [(pos(3) - SCATTER_LOGLOG_XTRIM) pos(4)]);
    set(h_fig, 'PaperPositionMode', 'manual');
    set(h_fig, 'PaperPosition',[(0 - .25*SCATTER_LOGLOG_XTRIM) 0 pos(3) pos(4)]);
    % set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
    print(h_fig, ...
          fullfile(resultDir, ...
                   'scatter_all_loglog.pdf'), ...
          '-dpdf');
    
    %scatter pgffa vs utppgffa
    h_fig = figure('Visible', 'off');
    hold on
    set(gca, 'FontSize', SCATTER_AXES_SIZE)
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    
    h_utppgffa = scatter(x(:), meanRuntimeUTPPGFFA(:), SCATTER_MARKER_SIZE, ...
                         'o', ...
                         'LineWidth', SCATTER_LINE_WIDTH);
    
    h_pgffa = scatter(x(:), meanRuntimePGFFA(:), SCATTER_MARKER_SIZE, ...
                      'o', ...
                      'LineWidth', SCATTER_LINE_WIDTH);

    title('\Lambda\rho vs Runtime', 'FontSize', SCATTER_TITLE_FONT)
    legend([h_pgffa, h_utppgffa], {'PGFFA', 'UTPPGFFA'}, 'Location', 'northwest', 'FontSize', SCATTER_LEGEND_FONT)
    % legend('show')
    xlabel('\Lambda\rho', 'FontSize', SCATTER_LABEL_FONT)
    ylabel('Mean runtime (s)', 'FontSize', SCATTER_LABEL_FONT)

    set(h_fig, 'PaperUnits','centimeters');
    set(h_fig, 'Units','centimeters');
    pos=get(h_fig,'Position');
    set(h_fig, 'PaperSize', [(pos(3) - SCATTER_LOGLOG_XTRIM) pos(4)]);
    set(h_fig, 'PaperPositionMode', 'manual');
    set(h_fig, 'PaperPosition',[(0 - .25*SCATTER_LOGLOG_XTRIM) 0 pos(3) pos(4)]);
    % set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
    print(h_fig, ...
          fullfile(resultDir, ...
                   'scatter_pgf_loglog.pdf'), ...
          '-dpdf');
    
    %scatter pgffa only
    h_fig = figure('Visible', 'off');
    hold on
    set(gca, 'FontSize', SCATTER_AXES_SIZE)
    
    colOrd = get(gca, 'ColorOrder');
    
    h_pgffa = scatter(x(:), meanRuntimePGFFA(:), SCATTER_MARKER_SIZE, ...
                      'o', ...
                      'LineWidth', SCATTER_LINE_WIDTH, ...
                      'MarkerFaceColor', colOrd(2,:));

    title('\Lambda\rho vs Runtime', 'FontSize', SCATTER_TITLE_FONT)
    legend([h_pgffa], {'PGFFA'}, 'Location', 'northwest', 'FontSize', SCATTER_LEGEND_FONT)
    % legend('show')
    xlabel('\Lambda\rho', 'FontSize', SCATTER_LABEL_FONT)
    ylabel('Mean runtime (s)', 'FontSize', SCATTER_LABEL_FONT)

    set(h_fig, 'PaperUnits','centimeters');
    set(h_fig, 'Units','centimeters');
    pos=get(h_fig,'Position');
    set(h_fig, 'PaperSize', [(pos(3) - SCATTER_LINEAR_XTRIM) pos(4)]);
    set(h_fig, 'PaperPositionMode', 'manual');
    set(h_fig, 'PaperPosition',[(0 - .25*SCATTER_LINEAR_XTRIM) 0 (pos(3)-SCATTER_LINEAR_XTRIM/2) pos(4)]);
    % set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
    print(h_fig, ...
          fullfile(resultDir, ...
                   'scatter_pgffa.pdf'), ...
          '-dpdf');
    
    %scatter utppgffa only
    h_fig = figure('Visible', 'off');
    hold on
    set(gca, 'FontSize', SCATTER_AXES_SIZE)
    
    h_utppgffa = scatter(x(:), meanRuntimeUTPPGFFA(:), SCATTER_MARKER_SIZE, ...
                      'o', ...
                      'LineWidth', SCATTER_LINE_WIDTH);

    title('\Lambda\rho vs Runtime', 'FontSize', SCATTER_TITLE_FONT)
    legend([h_utppgffa], {'UTPPGFFA'}, 'Location', 'northwest', 'FontSize', SCATTER_LEGEND_FONT)
    % legend('show')
    xlabel('\Lambda\rho', 'FontSize', SCATTER_LABEL_FONT)
    ylabel('Mean runtime (s)', 'FontSize', SCATTER_LABEL_FONT)

    set(h_fig, 'PaperUnits','centimeters');
    set(h_fig, 'Units','centimeters');
    pos=get(h_fig,'Position');
    set(h_fig, 'PaperSize', [(pos(3) - SCATTER_LINEAR_XTRIM) pos(4)]);
    set(h_fig, 'PaperPositionMode', 'manual');
    set(h_fig, 'PaperPosition',[(0 - .25*SCATTER_LINEAR_XTRIM) 0 (pos(3)-SCATTER_LINEAR_XTRIM/2) pos(4)]);
    % set(h_fig, 'PaperPosition',[(0 + xtrim/2) 0 (pos(3) - xtrim/2) pos(4)]);
    print(h_fig, ...
          fullfile(resultDir, ...
                   'scatter_utppgffa.pdf'), ...
          '-dpdf');
end