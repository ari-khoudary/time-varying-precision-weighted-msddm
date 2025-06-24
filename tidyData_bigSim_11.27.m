%% load in .mat files
if ~exist('outDir', 'var')
    dataDir = 'sim_files/bigSim_11.27';
else
    dataDir = outDir;
end
files = (dir([pwd filesep dataDir filesep '*.mat']));
allData = load([dataDir filesep files(1).name]).data;
allData = repmat(allData, length(files), 1);
for i=2:length(files)
    allData(i) = load([dataDir filesep files(i).name]).data;
end

if isfield(allData, 'dvSlopes')==0
    computeSlopes;
end

%% make long structure (each trial is one row) -- this should be a function
longData = repmat(struct('cue', [], 'coherence', [], 'congruent', [], 'threshold', [], 'memoryThinning', [], ...
    'noise1Frames', [], 'signal1Onsets', [], 'signal1Frames', [], 'noise2Onsets', [], 'noise2Frames', [], 'signal2Onsets', [], 'signal2Frames', [], ...
    'rawChoice', [], 'forcedChoice', [], 'RT', [], 'startPoint1', [], 'startPoint2', [], 'dvSlope1', [], 'dvSlope2', [], 'dvSlope3', [], 'dvSlope4', [], 'totalSlope', [], 'totalSlope_thresholded', [], ...
    'subID', [], 'nTrial', [], 'nFrames', []), sum([allData.nTrial]), 1);
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

    sig1o = num2cell(allData(i).signal1Onsets);
    [longData(trialCounter(i)+1:trialCounter(i+1)).signal1Onsets] = sig1o{:};

    sig1 = num2cell(allData(i).signal1Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).signal1Frames] = sig1{:};

    noise2o = num2cell(allData(i).noise2Onsets);
    [longData(trialCounter(i)+1:trialCounter(i+1)).noise2Onsets] = noise2o{:};

    noise2 = num2cell(allData(i).noise2Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).noise2Frames] = noise2{:};

    sig2o = num2cell(allData(i).signal2Onsets);
    [longData(trialCounter(i)+1:trialCounter(i+1)).signal2Onsets] = sig2o{:};

    sig2 = num2cell(allData(i).signal2Frames);
    [longData(trialCounter(i)+1:trialCounter(i+1)).signal2Frames] = sig2{:};

    rawC = num2cell(allData(i).choices(:,1));
    [longData(trialCounter(i)+1:trialCounter(i+1)).rawChoice] = rawC{:};

    forcedC = num2cell(allData(i).choices(:,2));
    [longData(trialCounter(i)+1:trialCounter(i+1)).forcedChoice] = forcedC{:};

    rt = num2cell(allData(i).RT);
    [longData(trialCounter(i)+1:trialCounter(i+1)).RT] = rt{:};

    sp1 = num2cell(allData(i).startPoints(:,1));
    [longData(trialCounter(i)+1:trialCounter(i+1)).startPoint1] = sp1{:}; % DV at start of trial

    sp2 = num2cell(allData(i).startPoints(:,2));
    [longData(trialCounter(i)+1:trialCounter(i+1)).startPoint2] = sp2{:}; % DV at onset of second signal period

    dv1 = num2cell(allData(i).dvSlopes(:,1)); 
    [longData(trialCounter(i)+1:trialCounter(i+1)).dvSlope1] = dv1{:}; % slope of DV during first signal period

    dv2 = num2cell(allData(i).dvSlopes(:,2));
    [longData(trialCounter(i)+1:trialCounter(i+1)).dvSlope2] = dv2{:}; % slope of DV during first noise period

    dv3 = num2cell(allData(i).dvSlopes(:,3));
    [longData(trialCounter(i)+1:trialCounter(i+1)).dvSlope3] = dv3{:}; % slope of DV during second signal period

    dv4 = num2cell(allData(i).dvSlopes(:,4));
    [longData(trialCounter(i)+1:trialCounter(i+1)).dvSlope4] = dv4{:}; % slope of DV during second noise period

    dvAll = num2cell(allData(i).dvSlopes(:,5));
    [longData(trialCounter(i)+1:trialCounter(i+1)).totalSlope] = dvAll{:}; % slope of dv from timepoint1:nFrames

    dvAll_thresh = num2cell(allData(i).dvSlopes(:,6));
    [longData(trialCounter(i)+1:trialCounter(i+1)).totalSlope_thresholded] = dvAll_thresh{:}; % slope of dv from timepoint 1 to RT

end

clear sid trials cues cohs threshs thins frames cong noise1 sig1o sig1 noise2o noise2 sig2o sig2 rawC forcedC rt sp1 sp2 dv1 dv2 dv3 dv4 dvAll dvAll_thresh trialCounter

%% convert to table 
dataTable = struct2table(longData);

%% compute relative RT
% as individual vectors
dataTable.noise1RT = dataTable.RT < dataTable.signal1Onsets;
dataTable.signal1RT = (dataTable.RT >= dataTable.signal1Onsets & dataTable.RT < dataTable.noise2Onsets);
dataTable.noise2RT = (dataTable.RT >= dataTable.noise2Onsets & dataTable.RT < dataTable.signal2Onsets);
dataTable.signal2RT = dataTable.RT >= dataTable.signal2Onsets;

% as a single categorical vector
dataTable.relativeRT = cell(height(dataTable),1);
dataTable.relativeRT(dataTable.noise1RT==1) = {'noise1'};
dataTable.relativeRT(dataTable.signal1RT==1) = {'signal1'};
dataTable.relativeRT(dataTable.noise2RT==1) = {'noise2'};
dataTable.relativeRT(dataTable.signal2RT==1) = {'signal2'};

%% create noise quartile variables
% right now this code is redundant because the noise durations come from
% the same generating distribution. but this is in case future me/us decide
% we want to use different PDFs for the different epochs

noise1quartiles = quantile(dataTable.noise1Frames, 4);
noise2quartiles = quantile(dataTable.noise2Frames, 4);
dataTable.noise1Quartile = zeros(height(dataTable), 1);
dataTable.noise2Quartile = dataTable.noise1Quartile;

for i=1:4
    if i==1
        dataTable.noise1Quartile(dataTable.noise1Frames <= noise1quartiles(i)) = i;
        dataTable.noise2Quartile(dataTable.noise2Frames <= noise2quartiles(i)) = i;
    elseif i<4
        dataTable.noise1Quartile(dataTable.noise1Frames <= noise1quartiles(i) & dataTable.noise1Frames > noise1quartiles(i-1)) = i;
        dataTable.noise2Quartile(dataTable.noise2Frames <= noise2quartiles(i) & dataTable.noise1Frames > noise1quartiles(i-1)) = i;
    else
        dataTable.noise1Quartile(dataTable.noise1Frames >= noise1quartiles(i)) = i;
        dataTable.noise2Quartile(dataTable.noise2Frames >= noise2quartiles(i)) = i;
    end
end

%% configure table for plotting
dataTable.congruent = categorical(dataTable.congruent, [1, 0], {'congruent', 'incongruent'});
dataTable.congruent(dataTable.cue==0.5) = 'neutral';
maxThresh = max(unique(dataTable.threshold));

%% make summary tables
% compute sem based on number of simulated subjects
nSub = allData(1).nSub;
f_sem = @(x)std(x, 'omitnan')/nSub;

% summary stats for probabilistic stimuli 
stimuliSummary = groupsummary(dataTable, ...
    ["congruent", "cue", "coherence", "threshold", "memoryThinning"], ... % grouping vars
     {'mean', f_sem, 'median'}, ... % operation
     {'noise1Frames', 'signal1Onsets', 'signal1Frames', 'noise2Onsets', 'noise2Frames', 'signal2Onsets', 'signal2Frames'}); % dependent vars

% summary stats for behavior
% NOTE: RT median in this table is not a good summary; its biased by the
% high threshold trials when DV doesn't hit threshold (and thus RT=nFrames)
% subsequent plots of median RT by condition compute the statistic directly
% in the plotting call (stat_summary)
summaryTable = groupsummary(dataTable, ...
    ["congruent", "cue", "coherence", "threshold", "memoryThinning"], ...
    {'mean', f_sem, 'median', 'mode'}, ...
    {'rawChoice', 'forcedChoice', 'RT', 'startPoint1', 'startPoint2', 'dvSlope1', 'dvSlope2', 'dvSlope3', 'dvSlope4', 'noise1RT', 'signal1RT', 'noise2RT', 'signal2RT'});

% add sem intervals to summary table
summaryTable.upperCI_rawChoice = summaryTable.mean_rawChoice + summaryTable.fun1_rawChoice;
summaryTable.lowerCI_rawChoice = summaryTable.mean_rawChoice - summaryTable.fun1_rawChoice;
summaryTable.upperCI_forcedChoice = summaryTable.mean_forcedChoice + summaryTable.fun1_forcedChoice;
summaryTable.lowerCI_forcedChoice = summaryTable.mean_forcedChoice - summaryTable.fun1_forcedChoice;
summaryTable.upperCI_RT = summaryTable.mean_RT + summaryTable.fun1_RT;
summaryTable.lowerCI_RT = summaryTable.mean_RT - summaryTable.fun1_RT;

% compute number of missing responses in raw choice
summaryTable.nanRT = groupsummary(dataTable, ["congruent", "cue", "coherence", "threshold", "memoryThinning"], 'nummissing', 'rawChoice').nummissing_rawChoice;

% reverse sign of start point for incongruent trials
dataTable.dummyStart1 = dataTable.startPoint1;
dataTable{dataTable.congruent=='incongruent', 'dummyStart1'} = -1*dataTable{dataTable.congruent=='incongruent', 'dummyStart1'};
summaryTable.dummyMedian_startPoint1 = summaryTable.median_startPoint1;
summaryTable{summaryTable.congruent=='incongruent', 'dummyMedian_startPoint1'} = -1*summaryTable{summaryTable.congruent=='incongruent', 'dummyMedian_startPoint1'};

% and of dv slope
dataTable.dummySlope1 = dataTable.dvSlope1;
dataTable.dummySlope2 = dataTable.dvSlope2;
dataTable.dummySlope3 = dataTable.dvSlope3;
dataTable.dummySlope4 = dataTable.dvSlope4;
% mean slope
summaryTable.dummyMean_dvSlope1 = summaryTable.mean_dvSlope1;
summaryTable.dummyMean_dvSlope2 = summaryTable.mean_dvSlope2;
summaryTable.dummyMean_dvSlope3 = summaryTable.mean_dvSlope3;
summaryTable.dummyMean_dvSlope4 = summaryTable.mean_dvSlope4;
% sem of slope
summaryTable.dummyFun1_dvSlope1 = summaryTable.fun1_dvSlope1;
summaryTable.dummyFun1_dvSlope2 = summaryTable.fun1_dvSlope2;
summaryTable.dummyFun1_dvSlope3 = summaryTable.fun1_dvSlope3;
summaryTable.dummyFun1_dvSlope4 = summaryTable.fun1_dvSlope4;
% median slope
summaryTable.dummyMedian_dvSlope1 = summaryTable.median_dvSlope1;
summaryTable.dummyMedian_dvSlope2 = summaryTable.median_dvSlope2;
summaryTable.dummyMedian_dvSlope3 = summaryTable.median_dvSlope3;
summaryTable.dummyMedian_dvSlope4 = summaryTable.median_dvSlope4;

dataNames = {'dummySlope1', 'dummySlope2', 'dummySlope3', 'dummySlope4'};
summaryNames = {'dummyMean_dvSlope1', 'dummyMean_dvSlope2', 'dummyMean_dvSlope3', 'dummyMean_dvSlope4'};
funNames = {'dummyFun1_dvSlope1', 'dummyFun1_dvSlope2', 'dummyFun1_dvSlope3', 'dummyFun1_dvSlope4'};
medNames = {'dummyMedian_dvSlope1', 'dummyMedian_dvSlope2', 'dummyMedian_dvSlope3', 'dummyMedian_dvSlope4'};

for i = 1:length(dataNames)
    dataTable{dataTable.congruent=='incongruent', dataNames{i}} = -1*dataTable{dataTable.congruent=='incongruent', dataNames{i}};
    summaryTable{summaryTable.congruent=='incongruent', summaryNames{i}} = -1*summaryTable{summaryTable.congruent=='incongruent', summaryNames{i}};
    summaryTable{summaryTable.congruent=='incongruent', funNames{i}} = -1*summaryTable{summaryTable.congruent=='incongruent', summaryNames{i}};
    summaryTable{summaryTable.congruent=='incongruent', medNames{i}} = -1*summaryTable{summaryTable.congruent=='incongruent', medNames{i}};
end

%% time-varying summary tables

% summary stats for behavior
% NOTE: RT median in this table is not a good summary; its biased by the
% high threshold trials when DV doesn't hit threshold (and thus RT=nFrames)
% subsequent plots of median RT by condition compute the statistic directly
% in the plotting call (stat_summary)
tv_summaryTable = groupsummary(dataTable, ...
    ["congruent", "cue", "coherence", "threshold", "memoryThinning", "noise1Quartile", "noise2Quartile"], ...
    {'mean', f_sem, 'median', 'mode'}, ...
    {'rawChoice', 'forcedChoice', 'RT', 'startPoint1', 'startPoint2', 'dvSlope1', 'dvSlope2', 'dvSlope3', 'dvSlope4', 'noise1RT', 'signal1RT', 'noise2RT', 'signal2RT'});

% add sem intervals to summary table
tv_summaryTable.upperCI_rawChoice = tv_summaryTable.mean_rawChoice + tv_summaryTable.fun1_rawChoice;
tv_summaryTable.lowerCI_rawChoice = tv_summaryTable.mean_rawChoice - tv_summaryTable.fun1_rawChoice;
tv_summaryTable.upperCI_forcedChoice = tv_summaryTable.mean_forcedChoice + tv_summaryTable.fun1_forcedChoice;
tv_summaryTable.lowerCI_forcedChoice = tv_summaryTable.mean_forcedChoice - tv_summaryTable.fun1_forcedChoice;
tv_summaryTable.upperCI_RT = tv_summaryTable.mean_RT + tv_summaryTable.fun1_RT;
tv_summaryTable.lowerCI_RT = tv_summaryTable.mean_RT - tv_summaryTable.fun1_RT;

% compute number of missing responses in raw choice
tv_summaryTable.nanRT = groupsummary(dataTable, ["congruent", "cue", "coherence", "threshold", "memoryThinning", 'noise1Quartile', 'noise2Quartile'], 'nummissing', 'rawChoice').nummissing_rawChoice;

% reverse sign of start point for incongruent trials
dataTable.dummyStart1 = dataTable.startPoint1;
dataTable{dataTable.congruent=='incongruent', 'dummyStart1'} = -1*dataTable{dataTable.congruent=='incongruent', 'dummyStart1'};
tv_summaryTable.dummyMedian_startPoint1 = tv_summaryTable.median_startPoint1;
tv_summaryTable{tv_summaryTable.congruent=='incongruent', 'dummyMedian_startPoint1'} = -1*tv_summaryTable{tv_summaryTable.congruent=='incongruent', 'dummyMedian_startPoint1'};

% mean slope
tv_summaryTable.dummyMean_dvSlope1 = tv_summaryTable.mean_dvSlope1;
tv_summaryTable.dummyMean_dvSlope2 = tv_summaryTable.mean_dvSlope2;
tv_summaryTable.dummyMean_dvSlope3 = tv_summaryTable.mean_dvSlope3;
tv_summaryTable.dummyMean_dvSlope4 = tv_summaryTable.mean_dvSlope4;
% sem of slope
tv_summaryTable.dummyFun1_dvSlope1 = tv_summaryTable.fun1_dvSlope1;
tv_summaryTable.dummyFun1_dvSlope2 = tv_summaryTable.fun1_dvSlope2;
tv_summaryTable.dummyFun1_dvSlope3 = tv_summaryTable.fun1_dvSlope3;
tv_summaryTable.dummyFun1_dvSlope4 = tv_summaryTable.fun1_dvSlope4;
% median slope
tv_summaryTable.dummyMedian_dvSlope1 = tv_summaryTable.median_dvSlope1;
tv_summaryTable.dummyMedian_dvSlope2 = tv_summaryTable.median_dvSlope2;
tv_summaryTable.dummyMedian_dvSlope3 = tv_summaryTable.median_dvSlope3;
tv_summaryTable.dummyMedian_dvSlope4 = tv_summaryTable.median_dvSlope4;

for i = 1:length(dataNames)
    tv_summaryTable{tv_summaryTable.congruent=='incongruent', summaryNames{i}} = -1*tv_summaryTable{tv_summaryTable.congruent=='incongruent', summaryNames{i}};
    tv_summaryTable{tv_summaryTable.congruent=='incongruent', funNames{i}} = -1*tv_summaryTable{tv_summaryTable.congruent=='incongruent', summaryNames{i}};
    tv_summaryTable{tv_summaryTable.congruent=='incongruent', medNames{i}} = -1*tv_summaryTable{tv_summaryTable.congruent=='incongruent', medNames{i}};
end