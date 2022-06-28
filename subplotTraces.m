
ssd = load('results/0.65cue_0.50coh_0cong_1subs_10trials_0.80antCoh_bayesUpdate.mat');


%%
subj = 1;
trial = 4;

addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];
figure; 

for i=1:trial
    subplot(trial/2, trial/2, i)
    hold on;
    plot(squeeze(ssd.memoryEvidence(subj,i,:)), 'LineWidth',2)
    plot([addX;squeeze(ssd.visualEvidence(subj, i, :))], 'LineWidth',2)
    plot([addX;squeeze(ssd.fullEvidence(subj,i,:))], 'LineWidth',2)
    plot([1,140],[0,0],'k')
    plot([ssd.nSampMemory,ssd.nSampMemory],[min(allEvidence)*1.2, max(allEvidence*1.2)],'k--')
    string2 = sprintf('subj=%i, trial=%i', subj, i);
    title(string2);
    xlabel('time (a.u.)')
    ylabel('evidence (a.u.)')
    makePretty(15)
end
string = sprintf('%s', ssd.outfile);
string = regexprep(string, '[\\\^\_]','\\$0');
sgtitle(string);
legend({'memory','visual','combined'},'location','southeast', 'FontSize',8)
hold off;
