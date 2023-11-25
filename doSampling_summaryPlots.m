%% doSampling_summaryPlots
% requires the gramm toolbox for matlab bc i'm a ggplotter to my core

%% load in .mat files
dataDir = outDir;
files = (dir([pwd filesep dataDir filesep '*.mat']));
allData = load([outDir filesep files(1).name]).data;
allData = repmat(allData, length(files), 1);
for i=2:length(files)
    allData(i) = load([outDir filesep files(i).name]).data;
end

%% make long structure (each trial is one row)
longData = repmat(struct('subID', [], 'nTrial', [], 'cue', [], 'coherence', [], 'congruent', [], 'threshold', [], 'memoryThinning', [], 'nFrames', [], ...
    'noise1Frames', [], 'noise2Frames', [], 'noise2Onset', [], 'signal1Frames', [], 'rawChoice', [], 'forcedChoice', [], 'RT', []), sum([allData.nTrial]), 1);
trialCounter = [0 cumsum([allData.nTrial])];

for i=1:length(allData)
    sid = repmat(num2cell(allData(i).subID), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).subID] = sid{:};

    trials = repmat(num2cell(allData(i).nTrial), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).nTrial] = trials{:};

    cues = repmat(num2cell(allData(i).cue), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).cue] = cues{:};

    cohs = repmat(num2cell(allData(i).coherence), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).coherence] = cohs{:};

    threshs = repmat(num2cell(allData(i).threshold), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).threshold] = threshs{:};

    thins = repmat(num2cell(allData(i).memoryThinning), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).memoryThinning] = thins{:};

    frames = repmat(num2cell(allData(i).nFrames), allData(i).nTrial, 1);
    [longData(trialCounter(i)+1:trialCounter(i+1)).nFrames] = frames{:};

    cong = num2cell(allData(i).congruent);
    [longData(trialCounter(i)+1:trialCounter(i+1)).congruent] = cong{:};

    noise1 = num2cell(allData(i).noise1Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).noise1Frames] = noise1{:};

    noise2 = num2cell(allData(i).noise2Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).noise2Frames] = noise2{:};

    noise2o = num2cell(allData(i).noise2Onset);
    [longData(trialCounter(i)+1:trialCounter(i+1)).noise2Onset] = noise2o{:};

    sig1 = num2cell(allData(i).signal1Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).signal1Frames] = sig1{:};

    rawC = num2cell(allData(i).choices(:,1));
    [longData(trialCounter(i)+1:trialCounter(i+1)).rawChoice] = rawC{:};

    forcedC = num2cell(allData(i).choices(:,2));
    [longData(trialCounter(i)+1:trialCounter(i+1)).forcedChoice] = forcedC{:};

    rt = num2cell(allData(i).RT);
    [longData(trialCounter(i)+1:trialCounter(i+1)).RT] = rt{:};

end

clear sid trials cues cohs threshs thins frames cong noise1 noise2 noise2o sig1 rawC forcedC rt

% and convert to table 
dataTable = struct2table(longData);

%% configure table for plotting
dataTable.congruent = categorical(dataTable.congruent, [1, 0], {'congruent', 'incongruent'});
dataTable.congruent(dataTable.cue==0.5) = 'neutral';

%% make summary tables
% compute sem based on number of simulated subjects
nSub = allData(1).nSub;
f_sem = @(x)std(x)/nSub;
semTable = grpstats(dataTable, ["congruent", "cue", "coherence", "threshold", "memoryThinning"], f_sem);
summaryTable = grpstats(dataTable, ["congruent", "cue", "coherence", "threshold", "memoryThinning"]);

% add sem intervals to summary table
summaryTable.upperCI_rawChoice = summaryTable.mean_rawChoice + semTable.Fun1_rawChoice;
summaryTable.lowerCI_rawChoice = summaryTable.mean_rawChoice - semTable.Fun1_rawChoice;
summaryTable.upperCI_forcedChoice = summaryTable.mean_forcedChoice + semTable.Fun1_forcedChoice;
summaryTable.lowerCI_forcedChoice = summaryTable.mean_forcedChoice - semTable.Fun1_forcedChoice;
summaryTable.upperCI_RT = summaryTable.mean_RT + semTable.Fun1_RT;
summaryTable.lowerCI_RT = summaryTable.mean_RT - semTable.Fun1_RT;

% compute number of missing responses in raw choice
summaryTable.nanRT = groupsummary(dataTable, ["congruent", "cue", "coherence", "threshold", "memoryThinning"], 'nummissing', 'rawChoice').nummissing_rawChoice;

%% plot summary effects of threshold on accuracy & RT (gramm method)

% raw choice accuracy
f(1,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_rawChoice, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_rawChoice, 'ymax', summaryTable.upperCI_rawChoice);
f(1,1).set_names('x', 'coherence', 'y', 'proportion correct (+/- sem)', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'threhsold');
f(1,1).set_title(['raw choice, n=' num2str(nSub)]);
f(1,1).facet_grid([], num2cell(num2str(summaryTable.threshold), 2));
f(1,1).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(1,1).geom_point();
f(1,1).geom_interval('geom', 'errorbar');
f(1,1).geom_line();

% forced choice accuracy
f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_forcedChoice, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_forcedChoice, 'ymax', summaryTable.upperCI_forcedChoice);
f(1,2).set_names('x', 'coherence', 'y', 'proportion correct (+/- sem)', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'threhsold');
f(1,2).set_title(['forced choice (sign of accumulator), n=' num2str(nSub)]);
f(1,2).facet_grid([], num2cell(num2str(summaryTable.threshold), 2));
f(1,2).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(1,2).geom_point();
f(1,2).geom_interval('geom', 'errorbar');
f(1,2).geom_line();


% proportion of not hitting threshold
f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.nanRT,...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(2,1).set_names('x', 'coherence', 'y', 'count', 'color', 'cue', 'linestyle', 'congruent', 'column', 'threshold');
f(2,1).set_title('number of trials with NaN choice');
f(2,1).facet_wrap(num2cell(num2str(summaryTable.threshold), 2));
f(2,1).geom_point();
f(2,1).geom_line();

% RT summary
f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_RT, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_RT, 'ymax', summaryTable.upperCI_RT);
f(2,2).set_names('x', 'coherence', 'y', 'mean RT', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'threhsold');
f(2,2).set_title(['mean +/- sem RT, n=' num2str(nSub)]);
f(2,2).facet_wrap(num2cell(num2str(summaryTable.threshold), 2));
f(2,2).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(2,2).geom_abline('slope', 0, 'intercept', max(dataTable.nFrames), 'style', '-');
f(2,2).geom_point();
f(2,2).geom_interval('geom', 'errorbar');
f(2,2).geom_line();

figure('Name', 'Effects of threshold on accuracy & RT', 'WindowState', 'maximized');
f.draw();









%% plot summary effects of threshold & coherence (matlab base plotting)
% accFC = figure;
% accFC.Name = 'threshold & accuracy (forced choice)';
% accRC = figure;
% accRC.Name = 'threshold & accuracy (raw choice)';
% 
% rtFC = figure;
% rtFC.Name = 'threshold & RT (forced choice)';
% rtRC = figure;
% rtRC.Name = 'threshold & RT (raw choice)';
% counter = 0;
% 
% for i = 1:length(coherenceLevels)
%     for j = 1:length(thresholdLevels)
%         titleString = sprintf('coherence=%.2f, threshold=%i', coherenceLevels(i), thresholdLevels(j));
%         for k = 1:length(cueLevels)
%             if cueLevels(k) == 0.5
%                 forcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j), :);
%                 choiceSEM = std(forcedChoice) / sqrt(length(unique(dataTable.subID)));
%                 rt = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j), :);
%                 rtSEM = std(rt) / sqrt(length(unique(dataTable.subID)));
% 
%                 figure(accFC)
%                 subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
%                 hold on
%                 ylim([0 1.5])
%                 xlim([0.3 1])
%                 xticks([0.3:0.1:1])
%                 yticks([0:0.1:1.25])
%                 yline(0.5, '--')
%                 yline(1, '-')
%                 grid on
%                 errorbar(cueLevels(k), mean(forcedChoice), choiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'black')
%                 xlabel('cue')
%                 ylabel('proportion correct')
%                 title(titleString);
% 
%                 figure(rtFC)
%                 subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
%                 hold on
%                 ylim([min(rt)-20 max(rt)+20])
%                 xlim([0.3 1])
%                 xticks([0.3:0.1:1])
%                 grid on
%                 errorbar(cueLevels(k), mean(rt), rtSEM, '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'black')
%                 xlabel('cue')
%                 ylabel('RT')
%                 title(titleString)
%             else
%                 congForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==1, :);
%                 congRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==1, :);
%                 incongForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==0, :);
%                 incongRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==0, :);
%                 congChoiceSEM = std(congForcedChoice) / sqrt(length(unique(dataTable.subID)));
%                 incongChoiceSEM = std(incongForcedChoice) / sqrt(length(unique(dataTable.subID)));
%                 congRTSEM = std(congRT) / sqrt(length(unique(dataTable.subID)));
%                 incongRTSEM = std(incongRT) / sqrt(length(unique(dataTable.subID)));
% 
%                 figure(accFC)
%                 h=subplot(length(thresholdLevels), length(coherenceLevels), counter+1);
%                 hold on
%                 errorbar(cueLevels(k), mean(congForcedChoice), congChoiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', "#7E2F8E")
%                 errorbar(cueLevels(k), mean(incongForcedChoice), incongChoiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', "#77AC30")
%                 if i==1 && j==1
%                     legend({'', '', 'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest');
%                 end
% 
%                 figure(rtFC)
%                 subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
%                 hold on
%                 errorbar(cueLevels(k), mean(congRT), congRTSEM, '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', "#7E2F8E")
%                 errorbar(cueLevels(k), mean(incongRT), incongRTSEM, 'x', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 	"#77AC30")
%                 if i==1 && j==1
%                     legend({'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest')
%                 end
% 
%             end
%         end
%         counter=counter+1;
%     end
% end
% 