function data = doSampling(config)

%% unpack params
nTrial = config.nTrial; %per cue
trialDuration = config.trialDuration; %in seconds; maximum possible duration
cue = config.cue;
coherence = config.coherence;
threshold = config.threshold;
memoryThinning = config.memoryThinning;
visionThinning = config.visionThinning;
vizPresentationRate = config.vizPresentationRate;
maxNoiseDuration = config.maxNoiseDuration;
minNoiseDuration = config.minNoiseDuration;
minSignalDuration = config.minSignalDuration;
flickerNoisePadding = config.flickerNoisePadding;
flickerNoiseValue = config.flickerNoiseValue;

% do you want to save any frame-by-frame information?
saveEvidence = config.saveEvidence;
saveAccumulators = config.saveAccumulators;
saveCounters = config.saveCounters;
savePrecisions = config.savePrecisions;
saveDrifts = config.saveDrifts;

% translate parameters into simulation-space
nFrames = trialDuration/vizPresentationRate;
congruentTrials = ceil(nTrial*cue);
congruent = [ones(congruentTrials, 1); zeros(nTrial-congruentTrials,1)];
noiseFrames = zeros(nTrial, 1);
maxNoiseFrames = maxNoiseDuration/vizPresentationRate;
noiseMin = round(minNoiseDuration/vizPresentationRate); %in frames
signalMin = round(minSignalDuration/vizPresentationRate); %in frames

%% create arrays to hold values
choices = zeros(nTrial, 2);
RTs = zeros(nTrial, 1);
memoryAccumulator = zeros(nFrames, nTrial);
visionAccumulator = zeros(nFrames, nTrial);
decisionVariable = zeros(nFrames, nTrial);
memoryPrecisions = zeros(nFrames, nTrial);
visionPrecisions = zeros(nFrames, nTrial);
memoryDrift = zeros(nFrames, nTrial);
visionDrift = zeros(nFrames, nTrial);
counters = zeros(nFrames, 4, nTrial);

% create index variables
alphaMem_idx = 1;
betaMem_idx = 2;
alphaVis_idx = 3;
betaVis_idx = 4;

%% set up
% precompute entropy distribution
p = 0.01:.01:0.99;
entropy = zeros(length(p), 1);
for i = 1:length(p)
    entropy(i) = computeEntropy(p(i));
end

% and exponential distribution for noise durations
lambda = 0.15;
noisePDF = discrete_bounded_hazard_rate(lambda, maxNoiseFrames);
noiseDistribution = round(noisePDF * nTrial);
trialCount = cumsum(noiseDistribution);
for i = 1:length(noiseDistribution)
    if i==1
        noiseFrames(1:noiseDistribution(i))=1;
    else
        startrow = trialCount(i) - noiseDistribution(i);
        endrow = trialCount(i);
        noiseFrames(startrow:endrow) = i;
    end
end

% create memory evidence stream
memoryStream = (binornd(1, cue, [nFrames, nTrial])*2-1) + normrnd(0,1, [nFrames, nTrial]);

% create visual evidence stream
noise1Frames = Shuffle(noiseFrames + noiseMin);
noise2Frames = Shuffle(noise1Frames);
signal1Frames = Shuffle(noiseFrames + signalMin);
noise2Onset = noise1Frames + signal1Frames;

if strcmp(flickerNoiseValue,'zero')==1
    flickerNoise = zeros(nFrames, nTrial);
elseif strcmp(flickerNoiseValue, 'gaussian')==1
    flickerNoise = normrnd(0,1, [nFrames, nTrial]);
end

if flickerNoisePadding==0
    flickerStream = ones(nFrames/2, nTrial);
else
    flickerStream = repmat([1; NaN], nFrames/2, nTrial);
    noiseBool = isnan(flickerStream);
    flickerStream(noiseBool) = flickerNoise(noiseBool);
end

for trial=1:nTrial
    flickerStream(1:noise1Frames(trial), trial)= flickerNoise(1:noise1Frames(trial));
    flickerStream(noise2Onset(trial):noise2Onset(trial)+noise2Frames(trial), trial) = flickerNoise(noise2Onset(trial):noise2Onset(trial)+noise2Frames(trial));
    imgIdx = find(flickerStream(trial,:)==1);
    nImgFrames = length(imgIdx);
    nTargetFrames = ceil(nImgFrames*coherence);
    targetIdx = randsample(imgIdx, nTargetFrames);
    lureIdx = setdiff(imgIdx, targetIdx);
    if trial <= congruentTrials
        target=1;
    else
        target=-1;
    end
    flickerStream(targetIdx, trial) = target;
    flickerStream(lureIdx, trial) = -target;
end

% compute analytic solution for each trial
expectedCounters = zeros(nTrial, 4);
expectedPrecisions = zeros(nTrial, 2);

expectedCounters(:, alphaMem_idx) = repmat(ceil(cue * (nFrames/memoryThinning)), [nTrial, 1]);
expectedCounters(:, betaMem_idx) = repmat((nFrames/memoryThinning), [nTrial,1]) - expectedCounters(:, alphaMem_idx);
if flickerNoisePadding==1
    expectedCounters(:, alphaVis_idx) = ceil(coherence * (nFrames/2 - (noise1Frames+noise2Frames)));
    expectedCounters(:, betaVis_idx) = (nFrames/2 - (noise1Frames+noise2Frames)) - expectedCounters(:, alphaVis_idx);
else
    expectedCounters(:, alphaVis_idx) = ceil(coherence * (nFrames - (noise1Frames+noise2Frames)));
    expectedCounters(:, betaVis_idx) = (nFrames - (noise1Frames+noise2Frames)) - expectedCounters(:, alphaVis_idx);
end

expectedPrecisions(:, 1) = 1./(betaVar(expectedCounters(:, alphaMem_idx), expectedCounters(:, betaMem_idx)));
expectedPrecisions(:, 2) = 1./(betaVar(expectedCounters(:, alphaVis_idx), expectedCounters(:, betaVis_idx)));

expectedAccuracy = expectedPrecisions(:, 2)>expectedPrecisions(:,1); 

%% run simulation

for trial=1:nTrial
    % initialize counters
    alphaMem = 1;
    betaMem = 1;
    alphaVis = 1;
    betaVis = 1;

    % compute time-varying drift rate
    for frame=1:nFrames

        if mod(frame, visionThinning)==0
            % update visual counters
            if flickerStream(frame, trial) > 0
                alphaVis = alphaVis + 1;
            elseif flickerStream(frame, trial) < 0
                betaVis = betaVis + 1;
            end

            % compute normalized beta pdf
            visualPDF = betapdf(p, alphaVis, betaVis) ./ sum(betapdf(p, alphaVis, betaVis));
            % compute precision as inverse belief-weighted entropy
            visualPrecision = 1/sum(visualPDF .* entropy');
            visionPrecisions(frame, trial) = visualPrecision;
        end

        % store counter values
        counters(frame, alphaVis_idx, trial) = alphaVis;
        counters(frame, betaVis_idx, trial) = betaVis;

        if mod(frame, memoryThinning)==1 % this allows a memory sample on the first frame
            % update memory counters
            if memoryStream(frame, trial) > 0
                alphaMem = alphaMem + 1;
            else
                betaMem = betaMem + 1;
            end

            % compute normalized beta pdf
            memoryPDF = betapdf(p, alphaMem, betaMem) ./ sum(betapdf(p, alphaMem, betaMem));
            % compute precision as inverse belief-weighted entropy
            memoryPrecision = 1/sum(memoryPDF .* entropy');
            memoryPrecisions(frame, trial) = memoryPrecision;
        end

        % store counter values
        counters(frame, alphaMem_idx, trial) = alphaMem;
        counters(frame, betaMem_idx, trial) = betaMem;

        % compute relative precision evidence weights
        visualDriftRate = visualPrecision / (visualPrecision + memoryPrecision);
        memoryDriftRate = memoryPrecision / (visualPrecision + memoryPrecision);
        visionDrift(frame, trial) = visualDriftRate;
        memoryDrift(frame, trial) = memoryDriftRate;

        % compute time-varying relative precision-weighted decision variable
        memorySample = memoryStream(frame, trial);
        visualSample = flickerStream(frame, trial);

        if frame == 1
            visionAccumulator(frame, trial) = visualSample * visualDriftRate;
            if mod(frame, memoryThinning)==1
                memoryAccumulator(frame, trial) = memorySample * memoryDriftRate;
                decisionVariable(frame, trial) = memorySample*memoryDriftRate + visualSample * visualDriftRate;
            else
                decisionVariable(frame, trial) = visualSample*visualDriftRate;
            end

        else % for frames > 1
            visionAccumulator(frame, trial) = visionAccumulator(frame-1, trial) + visualSample * visualDriftRate;
            if mod(frame, memoryThinning) == 1
                memoryAccumulator(frame, trial) = memoryAccumulator(frame-1, trial) + memorySample * memoryDriftRate;
                decisionVariable(frame, trial) = decisionVariable(frame-1, trial) + memorySample*memoryDriftRate + visualSample * visualDriftRate;
            else
                memoryAccumulator(frame, trial) = memoryAccumulator(frame-1, trial);
                decisionVariable(frame, trial) = decisionVariable(frame-1, trial) + visualSample*visualDriftRate;
            end
        end
    end

    % make decision
    % find point at which evidence crosses threshold
    boundaryIdx = find(abs(decisionVariable(:, trial))>threshold, 1);
    if isempty(boundaryIdx)
        boundaryIdx = nFrames;
    end

    % store crossing point as RT
    RTs(trial) = boundaryIdx;

    % populate choice matrix accordingly - first column is "raw"
    % responses, second is "forced"
    if trial <= congruentTrials
        if decisionVariable(boundaryIdx, trial) > threshold
            choices(trial, 1) = 1;
            choices(trial, 2) = 1;
        elseif boundaryIdx==nFrames
            choices(trial, 1) = NaN;
            if decisionVariable(boundaryIdx, trial) > 0
                choices(trial, 2) = 1;
            else
                choices(trial, 2) = 0;
            end
        else % if wrong bound is reached
            choices(trial, 1) = 0;
            choices(trial, 2) = 0;
        end
    else % if incongruent
        if decisionVariable(boundaryIdx, trial) < -threshold
            choices(trial, 1) = 1;
            choices(trial, 2) = 1;
        elseif boundaryIdx==nFrames
            choices(trial, 1) = NaN;
            if decisionVariable(boundaryIdx, trial) < 0
                choices(trial, 2) = 1;
            else
                choices(trial, 2) = 0;
            end
        else % if wrong bound is reached
            choices(trial, 1) = 0;
            choices(trial, 2) = 0;
        end
    end
end

%% store results
% simulation settings
data.nTrial = nTrial;
data.trialDuration = trialDuration;
data.cue = cue;
data.coherence = coherence;
data.threshold = threshold;
data.congruent = congruent;
data.memoryThinning = memoryThinning;
data.visionThinning = visionThinning;
data.nFrames = nFrames;
data.vizPresentationRate = vizPresentationRate;
data.maxNoiseDuration = maxNoiseDuration;
data.minNoiseDuration = minNoiseDuration;
data.minSignalDuration = minSignalDuration;
data.noise1Frames = noise1Frames;
data.noise2Frames = noise2Frames;
data.noise2Onset = noise2Onset;
data.signal1Frames = signal1Frames;
data.flickerNoiseValue = flickerNoiseValue;
data.flickerNoisePadding = flickerNoisePadding;
data.counters = counters;

% behavior
data.choices = choices;
data.RT = RTs;
data.expectedAccuracy = expectedAccuracy;

% optional frame-by-frame info
if saveEvidence==1
    data.memoryEvidence = memoryStream;
    data.visionEvidence = flickerStream;
end

if saveAccumulators==1
    data.memoryAccumulator = memoryAccumulator;
    data.visionAccumulator = visionAccumulator;
end

if saveCounters==1
    data.counters=counters;
end

if savePrecisions==1
    data.memoryPrecisions = memoryPrecisions;
    data.visionPrecisions = visionPrecisions;
end

if saveDrifts==1
    data.memoryDrifts = memoryDrift;
    data.visionDrifts = visionDrift;
end

end



