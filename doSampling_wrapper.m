%% specify simulation settings

nTrial = 10; % per cue
cue = [0.5];
coherence = [0.52];
threshold = [3];
memoryThinning = 12;
visionThinning = 1;
vizPresentationRate = 1/60;

% durations in seconds
expLambda = 0.15;
maxNoiseDuration = 2;
minNoiseDuration = 1;
minSignalDuration = 0.375;
secondSignalMin = 0; % value to be added to signalMin to create a second signal period of additional length

% half neutral trials (boolean)
halfNeutralTrials = 1;

% visual evidence noise
flickerNoisePadding = 1; %logical; do you want to pad each signal frame with a noise frame?
flickerNoiseValue = 'zero'; %string; how do you want to model noise frames? options are zeros, zero-centered gaussian, more to come

% do you want to save frame-by-frame information for each trial?
saveEvidence = 1;
saveFlickerNoise = 1;
saveAccumulators = 1;
saveCounters = 1;
savePrecisions = 1;
saveDrifts = 1;

%% create config file
config.nTrial = nTrial;
config.cue = cue;
config.coherence = coherence;
config.threshold = threshold;
config.memoryThinning = memoryThinning;
config.visionThinning = visionThinning;
config.vizPresentationRate = vizPresentationRate;

% durations in seconds
config.maxNoiseDuration = maxNoiseDuration;
config.minNoiseDuration = minNoiseDuration;
config.minSignalDuration = minSignalDuration;
config.secondSignalMin = secondSignalMin;
config.expLambda = expLambda;

config.halfNeutralTrials = halfNeutralTrials;
config.flickerNoisePadding = flickerNoisePadding;
config.flickerNoiseValue = flickerNoiseValue;
config.saveEvidence = saveEvidence;
config.saveFlickerNoise = saveFlickerNoise;
config.saveAccumulators = saveAccumulators;
config.saveCounters = saveCounters;
config.savePrecisions = savePrecisions;
config.saveDrifts = saveDrifts;

%% run simulation
tic
data = doSampling(config);
toc