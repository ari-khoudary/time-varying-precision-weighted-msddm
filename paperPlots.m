%% weak & strong cues with chance coherence
% weak first
ssd = load('results/0.65cue_0.50coh_0cong_100subs_25trials.mat');
subj = 1;
trial = 4;
addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];
fig=figure;

subplot(2,2,1)
plot(squeeze(ssd.memoryEvidence(subj,trial,:)), 'LineWidth',2)
hold on;
plot([addX;squeeze(ssd.visualEvidence(subj, trial, :))], 'LineWidth',2)
plot([addX;squeeze(ssd.fullEvidence(subj,trial,:))], 'LineWidth',2)
plot([1,140],[0,0],'k')
plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*1.2, max(allEvidence*1.2)],'k--')
title('weak cue, chance coherence')

% then strong
ssd = load('results/0.80cue_0.50coh_0cong_100subs_25trials.mat');
subj = 1;
trial = 10;
addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];

subplot(2,2,2)
plot(squeeze(ssd.memoryEvidence(subj,trial,:)), 'LineWidth',2)
hold on;
plot([addX;squeeze(ssd.visualEvidence(subj, trial, :))], 'LineWidth',2)
plot([addX;squeeze(ssd.fullEvidence(subj,trial,:))], 'LineWidth',2)
plot([1,140],[0,0],'k')
plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*1.2, max(allEvidence*1.2)],'k--')
title('strong cue, chance coherence')

% plots=axes(fig, 'visible', 'off');
% plots.XLabel.Visible='on';
% plots.YLabel.Visible='on';
% xlabel(plots, 'time (a.u.)');
% ylabel(plots, 'evidence (a.u.)');
% makePretty(15)


% congruent first
ssd = load('results/0.80cue_0.65coh_1cong_100subs_25trials.mat');
subj = 1;
trial = 4;
addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];

subplot(2,2,3)
plot(squeeze(ssd.memoryEvidence(subj,trial,:)), 'LineWidth',2)
hold on;
plot([addX;squeeze(ssd.visualEvidence(subj, trial, :))], 'LineWidth',2)
plot([addX;squeeze(ssd.fullEvidence(subj,trial,:))], 'LineWidth',2)
plot([1,140],[0,0],'k')
plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*5, max(allEvidence*1.4)],'k--')
subtitle('congruent')

% then incongruent
ssd = load('results/0.80cue_0.65coh_0cong_100subs_25trials.mat');
subj = 1;
trial = 6;
addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];

subplot(2,2,4)
plot(squeeze(ssd.memoryEvidence(subj,trial,:)), 'LineWidth',2)
hold on;
plot([addX;squeeze(ssd.visualEvidence(subj, trial, :))], 'LineWidth',2)
plot([addX;squeeze(ssd.fullEvidence(subj,trial,:))], 'LineWidth',2)
plot([1,140],[0,0],'k')
plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*2, max(allEvidence*1.5)],'k--')
subtitle('incongruent')
legend({'memory','visual','combined'},'location','northwest')

sgtitle('strong cue, weak coherence', 'FontWeight', 'bold');
plots=axes(fig, 'visible', 'off');
plots.XLabel.Visible='on';
plots.YLabel.Visible='on';
plots.Title.Visible='on';
xlabel(plots, 'time (a.u.)');
ylabel(plots, 'evidence (a.u.)');
makePretty(15)