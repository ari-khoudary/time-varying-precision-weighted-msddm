
ssd = load('results/0.80cue_0.80coh_0cong_1subs_10trials_0.80antCoh_bayesUpdate.mat');

%%
subj = 1;
trial = 4;

addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];

figure; hold on;
string = sprintf('%s', ssd.outfile);
string = regexprep(string, '[\\\^\_]','\\$0');
string2 = sprintf('subj=%i, trial=%i', subj, trial);
title(string);
subtitle(string2);
plot(squeeze(ssd.memoryEvidence(subj,trial,:)), 'LineWidth',2)
plot([addX;squeeze(ssd.visualEvidence(subj, trial, :))], 'LineWidth',2)
plot([addX;squeeze(ssd.fullEvidence(subj,trial,:))], 'LineWidth',2)
plot([1,140],[0,0],'k')
plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*1.2, max(allEvidence*1.2)],'k--')
xlabel('time (a.u.)')
ylabel('evidence (a.u.)')
legend({'memory','visual','combined'},'location','northwest')
makePretty(15)