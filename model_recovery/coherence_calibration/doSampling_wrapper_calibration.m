%% specify simulation settings
clear
nSub = 500; % per cell
nTrial = 1000; % per subject
coherence = 0.51:0.01:0.65;
threshold = 3:3:30;
memoryThinning = 12;
visionThinning = 1;
vizPresentationRate = 1/60;

% noise periods
noisePeriods = 0; % logical: do you want 2 noise periods?
noNoiseTrialDuration = 3; % in seconds: how long do you want trials to be if there are no noise periods?
% parameters of the noise periods
expLambda = 0.12; % parameter of the exponential defining the hazard rates 
maxNoiseDuration = 1.4; % seconds
minNoiseDuration = 0.75; % seconds
minSignalDuration = 0.5; % seconds
secondSignalMin = 0.5; % value to be added to signalMin to create a second signal period of additional length

% visual evidence noise
flickerAdditiveNoise = 1;  % logical; do you want to add noise to each sample of visual evidence?
flickerAdditiveNoiseValue = 'gaussian';  % string; what kind of noise do you want to add to each visual evidence sample? gaussian=zero-centered gaussian
flickerPadding = 1;  % logical; do you want to pad each signal frame with a noise frame?
flickerPaddingValue = 'zero'; % string; how do you want to model noise frames? options are zeros, zero-centered gaussian, more to come

% do you want to save frame-by-frame information for each trial?
saveEvidence = 0;
saveFlickerNoise = 0;
saveAccumulators = 0;
saveDV = 0;
saveCounters = 0;
savePrecisions = 0;
saveDrifts = 0;

% where do you want to save the results? (subdirectory of current dir)
outDir = 'results/';

%% create cell array to store config files
nCombo = length(coherence) * length(threshold);
allConfigs = repmat({struct('myfield', {})}, 1, nCombo);

% creates one config file for each combination of threshold & coherence 
counter=0;
for a = 1:length(coherence)
    for d = 1:length(threshold)

        counter=counter+1;

        config.nTrial = nTrial;
        config.coherence = coherence(a);
        config.threshold = threshold(d);
        config.vizPresentationRate = vizPresentationRate;
        config.noNoiseTrialDuration = noNoiseTrialDuration;

        config.flickerAdditiveNoise = flickerAdditiveNoise;
        config.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
        config.flickerNoisePadding = flickerPadding;
        config.flickerPaddingValue = flickerPaddingValue;
        config.saveEvidence = saveEvidence;
        config.saveFlickerNoise = saveFlickerNoise;
        config.saveAccumulators = saveAccumulators;
        config.saveDV = saveDV;
        config.saveCounters = saveCounters;
        config.savePrecisions = savePrecisions;
        config.saveDrifts = saveDrifts;
        config.outDir = outDir;

        allConfigs{counter} = config;
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
for c = 1:length(allConfigs)
    for subj = 1:nSub
        configIdx(counter) = c;
        subIdx(counter) = subj;
        counter = counter + 1;
    end
end

% start parallel pool
if isempty(gcp('nocreate'))
    parpool('local');
end

%% loop over iterations in parallel
parfor it = 1:totalIterations
    thisConfig = allConfigs{configIdx(it)};
    thisConfig.subID = subIdx(it);
    allData{it} = doSampling_calibration(thisConfig);
end



%% concatenate all tables
data = vertcat(allData{:});

% close parallel pool
delete(gcp('nocreate'));

%% write out data as a csv

if ~exist(outDir, 'dir')
    mkdir(outDir);
end
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
filename = sprintf('results_%s.csv', timestamp);
writetable(data, fullfile(outDir, filename));
