ssd = load('results/0.80cue_0.80coh_0cong_1subs_10trials_0.80antCoh_bayesUpdate.mat');

%%
subj = 1;
trial = 4;

addX = NaN(ssd.nSampMemory,1);
allEvidence = [squeeze(ssd.memoryEvidence(subj,trial,:));squeeze(ssd.visualEvidence(subj,trial,:));squeeze(ssd.fullEvidence(subj,trial,:))];
figure; 

for i=1:trial
    subplot(trial/2, trial/2, i)
    hold on;
    plot(squeeze(ssd.memoryDriftRates(subj,i,:)), 'LineWidth',2)
    plot([addX;squeeze(ssd.visualDriftRates(subj, i, :))], 'LineWidth',2)
    plot([1,140],[0,0],'k')
    xlabel('time (a.u.)')
    ylabel('drift rate (a.u.)')
    string2 = sprintf('subj=%i, trial=%i', subj, i);
    title(string2);
    makePretty(15)
end
legend({'memory','visual'},'location','southeast')
string = sprintf('%s', ssd.outfile);
string = regexprep(string, '[\\\^\_]','\\$0');
sgtitle(string);
hold off
