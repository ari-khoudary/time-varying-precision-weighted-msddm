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
%% add fields to be populated



