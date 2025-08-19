%% specify simulation settings
clear
debug = 0;
cue = [0.5 0.65 0.8];
threshold_coherence = readmatrix("calibratedCoherences.csv");
threshold_coherence = threshold_coherence(:, 1:2);
memoryThinning = 4:4:48;
visionThinning = 1;
vizPresentationRate = 1/60;
if ~debug
    nSub = 500;
    nTrial = 1000; 
else
    nSub = 1;
    nTrial = 10;
end

% noise periods -- correspond to expt design as of 08-2025
noisePeriods = 1; % logical: do you want 2 noise periods?
noNoiseTrialDuration = 3; % in seconds: how long do you want trials to be if there are no noise periods?
% parameters of the noise periods
expLambda = 0.12; % parameter of the exponential defining the hazard rates 
minNoiseDuration = 0.75; % seconds
minSignalDuration = 0.5; % seconds

% half neutral trials (boolean)
halfNeutralTrials = 1;

% visual evidence noise
flickerAdditiveNoise = 1;  % logical; do you want to add noise to each sample of visual evidence?
flickerAdditiveNoiseValue = 'gaussian';  % string; what kind of noise do you want to add to each visual evidence sample? gaussian=zero-centered gaussian
flickerPadding = 1;  % logical; do you want to pad each signal frame with a noise frame?
flickerPaddingValue = 'zero'; % string; how do you want to model noise frames? options are zeros, zero-centered gaussian, more to come

% where do you want to save the results? (subdirectory of current dir)
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd'));
outDir = sprintf('results/%s/', timestamp);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% create config files
nCombo = length(threshold_coherence) * length(memoryThinning) * length(cue);
allConfigs = repmat({struct('myfield', {})}, 1, nCombo);

counter=0;
for a = 1:length(threshold_coherence)
    for b = 1:length(cue)
        for s = 1:length(memoryThinning)

            counter=counter+1;

            config.nTrial = nTrial;
            config.threshold = threshold_coherence(a, 1);
            config.coherence = threshold_coherence(a, 2);
            config.cue = cue(b);
            config.memoryThinning = memoryThinning(s);
            config.visionThinning = visionThinning;
            config.vizPresentationRate = vizPresentationRate;
            config.noisePeriods = noisePeriods;
            config.noNoiseTrialDuration = noNoiseTrialDuration;

            %config.maxNoiseDuration = maxNoiseDuration;
            config.minNoiseDuration = minNoiseDuration;
            config.minSignalDuration = minSignalDuration;
            %config.secondSignalMin = secondSignalMin;
            config.expLambda = expLambda;

            config.halfNeutralTrials = halfNeutralTrials;
            config.flickerAdditiveNoise = flickerAdditiveNoise;
            config.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
            config.flickerNoisePadding = flickerPadding;
            config.flickerPaddingValue = flickerPaddingValue;
            config.outDir = outDir;

            allConfigs{counter} = config;
        end
    end
end

%% run simulation in parallel
% calculate total number of iterations
totalIterations = length(allConfigs) * nSub;

% initialize structure to store output of each iteraction
allData = cell(totalIterations, 1);

% create arrays for parallel processing
configIdx = zeros(totalIterations, 1);
subIdx = zeros(totalIterations, 1);

% populate indices
counter = 1;
for s = 1:length(allConfigs)
    for subj = 1:nSub
        configIdx(counter) = s;
        subIdx(counter) = subj;
        counter = counter + 1;
    end
end

% start parallel pool
if ~debug
    if isempty(gcp('nocreate'))
        parpool('local');
    end
end

%% loop over iterations in parallel

parfor it = 1:totalIterations
    subIteration = allConfigs{configIdx(it)};
    subIteration.subID = subIdx(it);
    allData{it} = doSampling(subIteration);
end

% save one csv per subject
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for s = 1:nSub
    subIts = find(subIdx == s);
    combinedData = [];

    for it = subIts'  % Loop through all iterations for this subject
        if isempty(combinedData)
            combinedData = allData{it};
        else
            combinedData = [combinedData; allData{it}];
        end
    end

    outfile = [outDir 'sub' num2str(s) '.csv'];
    writetable(combinedData, outfile);
end







