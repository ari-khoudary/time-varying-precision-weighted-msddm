% trial plots
% traces
if plotResults==1
    fig=figure;
    for a=1:10
        subplot(2, 5, a)
        hold on;
        plot(memoryAccumulator(:, a, 1), 'LineWidth',1.5)
        plot(visionAccumulator(:, a, 1), 'LineWidth',1.5)
        plot(decisionVariable(:, a, 1), 'LineWidth',1.5)
        xline(noise1Frames(a))
        xline([noise1Frames(a)+signal1Frames(a), noise1Frames(a)+signal1Frames(a)+noise2Frames(a)], 'k--')
        yline([threshold -threshold])
        string2 = sprintf('trial %i', a);
        title(string2);
    end

    % add title
    if congruent==1
        string = sprintf('%.2f cue, %.2f coherece, %i thinning, congruent', cueLevel, coherence, memoryThinning);
    else
        string = sprintf('%.2f cue, %.2f coherence, %i thinning, incongruent', cueLevel, coherence, memoryThinning);
    end
    sgtitle(sprintf(['combo model\n' string]));

    % make pretty
    h=legend({'memory','visual','combined'},'FontSize',8, 'Orientation', 'horizontal');
    set(h, 'Position', [0.65 0.46 0.25 0.025]);
    plots=axes(fig, 'visible', 'off');
    plots.XLabel.Visible='on';
    plots.YLabel.Visible='on';
    plots.Title.Visible='on';
    xlabel(plots, 'time (a.u.)');
    ylabel(plots, 'evidence (a.u.)');
    if writeCSV==1
        figpath = sprintf('results_v2/trace_figs/%.2fcue_%icong_%.2fcoh_%ithin.png', cueLevel, congruent, coherence, memoryThinning);
        saveas(gcf, figpath);
    end
end










