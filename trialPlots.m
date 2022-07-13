%% specify configuration & load in file
dir = 'bayes_results/';
cueLevel = 0.80;
anticipatedCoherence = 0.50;
coherenceLevel = 0.65;
congruent = 0;
trials=1:10;

infile = sprintf([dir '%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat'], cueLevel, anticipatedCoherence, coherenceLevel, congruent);
load(infile);
addX = NaN(nSampMemory,1);

%% traces
fig=figure; 
counter=1;
for a=trials
    subplot(2, 5, counter)
    hold on;
    plot(squeeze(memoryEvidence(subj,a,:)), 'LineWidth',1.5)
    plot([addX;squeeze(visualEvidence(subj, a, :))], 'LineWidth',1.5)
    plot([addX;squeeze(fullEvidence(subj,a,:))], 'LineWidth',1.5)
    plot([1,140],[0,0],'k')
    xline(nSampMemory, 'k--')
    yline([threshold -threshold])
    string2 = sprintf('trial %i', a);
    title(string2);
    counter=counter+1;
end
if congruent==1
    string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, congruent', cueLevel, anticipatedCoherence, coherenceLevel);
else
    string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, incongruent', cueLevel, anticipatedCoherence, coherenceLevel);
end
if strcmp(dir, 'entropy_results/')
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
%saveas(gcf, [outfig '_traces.png']);
 
%% drifts
fig=figure; 
counter=1;
for b=trials
    subplot(2, 5, counter)
    hold on;
    plot(squeeze(memoryDriftRates(subj,b,:)), 'LineWidth',1.5)
    plot([addX;squeeze(visualDriftRates(subj, b, :))], 'LineWidth',1.5)
    plot([1,140],[0,0],'k')
    xline(nSampMemory, 'k--')
    string2 = sprintf('trial %i', b);
    title(string2);
    counter=counter+1;
end

h=legend({'memory','visual'},'Orientation','horizontal');
set(h, 'Position', [0.65 0.46 0.25 0.025]);

if strcmp(dir, 'entropy_results/')
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
%saveas(gcf, [outfig '_traces.png']);

%% precisions
fig=figure; 
counter=1;
for b=trials
    subplot(2, 5, counter)
    hold on;
    plot(squeeze(memoryPrecisions(subj,b,:)), 'LineWidth',1.5)
    plot([addX;squeeze(visualPrecisions(subj, b, :))], 'LineWidth',1.5)
    plot([1,140],[0,0],'k')
    xline(nSampMemory, 'k--')
    string2 = sprintf('trial %i', b);
    title(string2);
    counter=counter+1;
end

h=legend({'memory','visual'},'Orientation','horizontal');
set(h, 'Position', [0.65 0.46 0.25 0.025]);

if strcmp(dir, 'entropy_results/')
     sgtitle(sprintf(['first-order model\n' string]));
else
     sgtitle(sprintf(['second-order model\n' string]));
end
plots=axes(fig, 'visible', 'off');
plots.XLabel.Visible='on';
plots.YLabel.Visible='on';
plots.Title.Visible='on';
xlabel(plots, 'time (a.u.)');
ylabel(plots, 'precision (a.u.)');
%saveas(gcf, [outfig '_traces.png']);