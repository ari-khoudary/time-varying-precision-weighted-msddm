ssd = load('results/0.80cue_0.80coh_0cong_1subs_10trials_bayesUpdate.mat');

%%
subj = 1;
trial = 10;

addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];

figure; hold on;
string = sprintf('%s', ssd.outfile);
string = regexprep(string, '[\\\^\_]','\\$0');
string2 = sprintf('subj=%i, trial=%i', subj, trial);
title(string);
subtitle(string2);
plot(squeeze(ssd.memoryPrecisions(subj,trial,:)), 'LineWidth',1)
plot([addX;squeeze(ssd.visualPrecisions(subj, trial, :))], 'LineWidth',1)
plot([1,140],[0,0],'k')
%plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*1.2, max(allEvidence*1.2)],'k--')
xlabel('time (a.u.)')
ylabel('precision (a.u.)')
ylim([1.4 2])
legend({'memory','visual'},'location','northwest')
makePretty(15)