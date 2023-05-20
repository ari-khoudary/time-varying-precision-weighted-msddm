%% specify simulation settings
nSub = 3;
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

% where do you want to save the results? (subdirectory of current dir)
outDir = 'v3/test';

% create cell array to store config files
nCombo = length(coherence)*length(cue)*length(threshold)*length(memoryThinning);
configs = repmat({struct('myfield', {})}, 1, nCombo);


%% create config files

counter=0;
for a = 1:length(coherence)
    for b = 1:length(cue)
        for c = 1:length(threshold)
            for d = 1:length(memoryThinning)

                counter=counter+1;
                
                config.nTrial = nTrial;
                config.nSub = nSub;
                config.coherence = coherence(a);
                config.cue = cue(b);
                config.threshold = threshold(c);
                config.memoryThinning = memoryThinning(d);
                config.visionThinning = visionThinning;
                config.vizPresentationRate = vizPresentationRate;
                
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
                config.outDir = outDir;

                configs{counter} = config;
            end
        end
    end
end

%% run simulation
tic
counter=0;
%data = repmat({struct('myfield', {})}, 1, nSub);
for a = 1:length(coherence)
    for b = 1:length(cue)
        for c = 1:length(threshold)
            for d = 1:length(memoryThinning)
                counter=counter+1;
                config = configs{counter};
                for subj = 1:nSub
                    config.subID = subj;
                    doSampling(config);
                end
            end
        end
    end
end
toc