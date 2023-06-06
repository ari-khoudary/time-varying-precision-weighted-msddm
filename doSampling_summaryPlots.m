%% doSampling_summaryPlots

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

% and convery to table 
dataTable = struct2table(longData);

%% make plotting variables
nCells = length(allData) / allData(1).nSub;
cueLevels = unique(dataTable.cue);
cueCat = categorical(cueLevels);
coherenceLevels = unique(dataTable.coherence);
thresholdLevels = unique(dataTable.threshold);
thinningLevels = unique(dataTable.memoryThinning);

%% plot summary effects of threshold & coherence
accPlot = figure;
accPlot.Name = 'threshold & accuracy';
rtPlot = figure;
rtPlot.Name = 'threshold & RT';
counter = 0;

for i = 1:length(coherenceLevels)
    for j = 1:length(thresholdLevels)
        titleString = sprintf('coherence=%.2f, threshold=%i', coherenceLevels(i), thresholdLevels(j));
        for k = 1:length(cueLevels)
            if cueLevels(k) == 0.5
                forcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j), :);
                choiceSEM = std(forcedChoice) / sqrt(length(unique(dataTable.subID)));
                rt = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j), :);
                rtSEM = std(rt) / sqrt(length(unique(dataTable.subID)));

                figure(accPlot)
                subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
                hold on
                ylim([0 1.5])
                xlim([0.3 1])
                xticks([0.3:0.1:1])
                yline(0.5, '--')
                yline(1, '-')
                grid on
                errorbar(cueLevels(k), mean(forcedChoice), choiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                xlabel('cue')
                ylabel('proportion correct')
                title(titleString);

                figure(rtPlot)
                subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
                hold on
                ylim([min(rt)-20 max(rt)+20])
                xlim([0.3 1])
                xticks([0.3:0.1:1])
                grid on
                errorbar(cueLevels(k), mean(rt), rtSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                xlabel('cue')
                ylabel('RT')
                title(titleString)
            else
                congForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==1, :);
                congRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==1, :);
                incongForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==0, :);
                incongRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.threshold==thresholdLevels(j) & dataTable.congruent==0, :);
                congChoiceSEM = std(congForcedChoice) / sqrt(length(unique(dataTable.subID)));
                incongChoiceSEM = std(incongForcedChoice) / sqrt(length(unique(dataTable.subID)));
                congRTSEM = std(congRT) / sqrt(length(unique(dataTable.subID)));
                incongRTSEM = std(incongRT) / sqrt(length(unique(dataTable.subID)));

                figure(accPlot)
                h=subplot(length(thresholdLevels), length(coherenceLevels), counter+1);
                hold on
                errorbar(cueLevels(k), mean(congForcedChoice), congChoiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                errorbar(cueLevels(k), mean(incongForcedChoice), incongChoiceSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                if i==1 && j==1
                    legend({'', '', 'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest');
                end

                figure(rtPlot)
                subplot(length(thresholdLevels), length(coherenceLevels), counter+1)
                hold on
                errorbar(cueLevels(k), mean(congRT), congRTSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                errorbar(cueLevels(k), mean(incongRT), incongRTSEM, '.', 'MarkerSize', 20, 'LineWidth', 2)
                if i==1 && j==1
                    legend({'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest')
                end

            end
        end
        counter=counter+1;
    end
end

%% plot summary effects of thinning & coherence
accPlot = figure;
accPlot.Name = 'thinning & accuracy';
rtPlot = figure;
rtPlot.Name = 'thinning & RT';
counter = 0;

for i = 1:length(coherenceLevels)
    for j = 1:length(thinningLevels)
        titleString = sprintf('coherence=%.2f, thinning=%i', coherenceLevels(i), thinningLevels(j));
        for k = 1:length(cueLevels)
            if cueLevels(k) == 0.5
                forcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j), :);
                choiceSEM = std(forcedChoice) / sqrt(length(unique(dataTable.subID)));
                rt = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j), :);
                rtSEM = std(rt) / sqrt(length(unique(dataTable.subID)));

                figure(accPlot)
                subplot(length(thinningLevels), length(coherenceLevels), counter+1)
                hold on
                ylim([0 1.5])
                xlim([0.3 1])
                xticks([0.3:0.1:1])
                yticks([0 0.25 0.5 0.75 1 1.25])
                yline(0.5, '--')
                yline(1, '-')
                grid on
                errorbar(cueLevels(k), mean(forcedChoice), choiceSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                xlabel('cue')
                ylabel('proportion correct')
                title(titleString);

                figure(rtPlot)
                subplot(length(thinningLevels), length(coherenceLevels), counter+1)
                hold on
                ylim([min(rt)-20 max(rt)+20])
                xlim([0.3 1])
                xticks([0.3:0.1:1])
                grid on
                errorbar(cueLevels(k), mean(rt), rtSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                xlabel('cue')
                ylabel('RT')
                title(titleString)
            else
                congForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j) & dataTable.congruent==1, :);
                congRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j) & dataTable.congruent==1, :);
                incongForcedChoice = dataTable.forcedChoice(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j) & dataTable.congruent==0, :);
                incongRT = dataTable.RT(dataTable.cue==cueLevels(k) & dataTable.coherence==coherenceLevels(i) & dataTable.memoryThinning==thinningLevels(j) & dataTable.congruent==0, :);
                congChoiceSEM = std(congForcedChoice) / sqrt(length(unique(dataTable.subID)));
                incongChoiceSEM = std(incongForcedChoice) / sqrt(length(unique(dataTable.subID)));
                congRTSEM = std(congRT) / sqrt(length(unique(dataTable.subID)));
                incongRTSEM = std(incongRT) / sqrt(length(unique(dataTable.subID)));

                figure(accPlot)
                h=subplot(length(thinningLevels), length(coherenceLevels), counter+1);
                hold on
                errorbar(cueLevels(k), mean(congForcedChoice), congChoiceSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                errorbar(cueLevels(k), mean(incongForcedChoice), incongChoiceSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                if i==1 && j==1
                    legend({'', '', 'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest');
                end

                figure(rtPlot)
                subplot(length(thinningLevels), length(coherenceLevels), counter+1)
                hold on
                errorbar(cueLevels(k), mean(congRT), congRTSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                errorbar(cueLevels(k), mean(incongRT), incongRTSEM, '*', 'MarkerSize', 10, 'LineWidth', 2)
                if i==1 && j==1
                    legend({'neutral', 'congruent', 'incongruent'}, 'Location', 'NorthWest')
                end

            end
        end
        counter=counter+1;
    end
end

%% plot dynamic effects of threshold & coherence

%% plot dynamic effects of thinning & coherence
