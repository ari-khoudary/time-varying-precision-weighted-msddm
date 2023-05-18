% plot traces & drifts 

dir = {'entropy_results/', 'bayes_results/'};
cues = [0.50, 0.65, 0.8];
antCohs = [0.50, 0.65, 0.8];
coherences = [0.50, 0.65, 0.8];
congruence = [0, 1];

for i = 1:length(dir)
    rootdir = dir{i};
    for j=1:length(cues)
        for k = 1:length(antCohs)
            for l = 1:length(coherences)
                for m = 1:length(congruence)
                    % load in the file
                    infile=sprintf('%s%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', rootdir, cues(j), antCohs(k), coherences(l), congruence(m));
                    load(infile);
    
                    % make plots
                    addX = NaN(nSampMemory,1);
                    fig=figure; 
                    for a=1:10
                        subplot(2, 5, a)
                        hold on;
                        plot(squeeze(memoryEvidence(subj,a,:)), 'LineWidth',1.5)
                        plot([addX;squeeze(visualEvidence(subj, a, :))], 'LineWidth',1.5)
                        plot([addX;squeeze(fullEvidence(subj,a,:))], 'LineWidth',1.5)
                        plot([1,140],[0,0],'k')
                        xline(nSampMemory, 'k--')
                        yline([threshold -threshold])
                        string2 = sprintf('trial %i', a);
                        title(string2);
                    end
                    if congruent
                        string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, congruent', cueLevel, anticipatedCoherence, coherenceLevel);
                    else
                        string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, incongruent', cueLevel, anticipatedCoherence, coherenceLevel);
                    end
    
                    if strcmp(rootdir, 'entropy_results/')
                         sgtitle(sprintf(['first-order model\n' string]));
                    else
                         sgtitle(sprintf(['second-order model\n' string]));
                    end
    
                    h=legend({'memory','visual','combined'},'FontSize',8, 'Orientation', 'horizontal');
                    set(h, 'Position', [0.55 0.46 0.35 0.025]);
                    outfig = char(regexp(outfile, '0.*cong', 'match'));
                    plots=axes(fig, 'visible', 'off');
                    plots.XLabel.Visible='on';
                    plots.YLabel.Visible='on';
                    plots.Title.Visible='on';
                    xlabel(plots, 'time (a.u.)');
                    ylabel(plots, 'evidence (a.u.)');
                    saveas(gcf, [rootdir 'plots/' outfig '_traces.png']);
                     
                    % drifts
                    fig=figure; 
                    for b=1:10
                        subplot(2, 5, b)
                        hold on;
                        plot(squeeze(memoryDriftRates(subj,b,:)), 'LineWidth',1.5)
                        plot([addX;squeeze(visualDriftRates(subj, b, :))], 'LineWidth',1.5)
                        plot([1,140],[0,0],'k')
                        xline(nSampMemory, 'k--')
                        string2 = sprintf('trial %i', b);
                        title(string2);
                    end
                    h=legend({'memory','visual'},'Orientation','horizontal');
                    set(h, 'Position', [0.65 0.46 0.25 0.025]);
    
                    if strcmp(rootdir, 'entropy_results/')
                         sgtitle(sprintf(['first-order model\n' string]));
                    else
                         sgtitle(sprintf(['second-order model\n' string]));
                    end
    
                    plots=axes(fig, 'visible', 'off');
                    plots.XLabel.Visible='on';
                    plots.YLabel.Visible='on';
                    plots.Title.Visible='on';
                    xlabel(plots, 'time (a.u.)');
                    ylabel(plots, 'drift rate (a.u.)');
                    saveas(gcf, [rootdir 'plots/' outfig '_drifts.png']);
    
                    close all
                end
            end
        end
    end
end


