%% specify simulation settings
clear
nSub = 1;
nTrial = 10; % per cue
cue = 0.8; % defines theta_memory
coherence = 0.5; % defines theta_vision
threshold = 15; % defines a in eq 7; irrelevant for all reported data
memoryThinning = [28]; % defines gamma in eq8, we tested at these levels
visionThinning = 1;
vizPresentationRate = 1/60; 
trialDuration = 3;

% early effect settings 
delayPeriod = 0; % logical: do you want a delay period between cue and visual evidence?
delayDurations = 0.75 + [4 6 8]; % duration of cue + possible ISI values to be drawn at uniform

% half neutral trials (boolean); only relevant if cue==0.5
halfNeutralTrials = 0;

% do you want to save frame-by-frame information for each trial?
saveEvidence = 1;
saveAccumulators = 1;
saveDV = 1;
saveCounters = 1;
savePrecisions = 1;
saveDrifts = 1;

% where do you want to save the results? (subdirectory of current dir)
outDir = 'extended_results';

%% create cell array to store config files
nCombo = length(coherence)*length(cue)*length(threshold)*length(memoryThinning);
allConfigs = repmat({struct('myfield', {})}, 1, nCombo);

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
                config.trialDuration = trialDuration;
                config.delayPeriod = delayPeriod;
                config.delayDurations = delayDurations;
                config.halfNeutralTrials = halfNeutralTrials;
                
                config.saveEvidence = saveEvidence;
                config.saveAccumulators = saveAccumulators;
                config.saveDV = saveDV;
                config.saveCounters = saveCounters;
                config.savePrecisions = savePrecisions;
                config.saveDrifts = saveDrifts;
                config.outDir = outDir;

                allConfigs{counter} = config;
            end
        end
    end
end

%% run simulation
tic
counter=0;
for a = 1:length(coherence)
    for b = 1:length(cue)
        for c = 1:length(threshold)
            for d = 1:length(memoryThinning)
                counter=counter+1;
                thisConfig = allConfigs{counter};
                sprintf(['running coherence=' num2str(coherence(a)) ', cue=' num2str(cue(b)) ', threshold=' num2str(threshold(c)) ' and thinning=' num2str(memoryThinning(d))])
                for subj = 1:nSub
                    thisConfig.subID = subj;
                    doSampling_ejn(thisConfig);
                end
            end
        end
    end
end
toc









