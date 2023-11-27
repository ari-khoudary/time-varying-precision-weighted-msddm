function data = doSampling(config)

%% unpack params
nTrial = config.nTrial; %per cue
nSub = config.nSub;
subID = config.subID;
cue = config.cue;
coherence = config.coherence;
threshold = config.threshold;
memoryThinning = config.memoryThinning;
visionThinning = config.visionThinning;
vizPresentationRate = config.vizPresentationRate;
noisePeriods = config.noisePeriods;
expLambda = config.expLambda;
maxNoiseDuration = config.maxNoiseDuration;
minNoiseDuration = config.minNoiseDuration;
minSignalDuration = config.minSignalDuration;
secondSignalMin = config.secondSignalMin;
halfNeutralTrials = config.halfNeutralTrials;
flickerAdditiveNoise = config.flickerAdditiveNoise;
flickerAdditiveNoiseValue = config.flickerAdditiveNoiseValue;
flickerNoisePadding = config.flickerNoisePadding;
flickerPaddingValue = config.flickerPaddingValue;

% do you want to save any frame-by-frame information?
saveEvidence = config.saveEvidence;
saveFlickerNoise = config.saveFlickerNoise;
saveAccumulators = config.saveAccumulators;
saveDV = config.saveDV;
saveCounters = config.saveCounters;
savePrecisions = config.savePrecisions;
saveDrifts = config.saveDrifts;

% where do you want to save the results?
outDir = config.outDir;

%% translate parameters into simulation-space
if cue==0.5 && halfNeutralTrials==1
    nTrial=nTrial/2;
end
congruentTrials = ceil(nTrial*cue);
congruent = [ones(congruentTrials, 1); zeros(nTrial-congruentTrials,1)];

% make noise durations
noiseFrames = zeros(nTrial, 1);
maxNoiseFrames = maxNoiseDuration/vizPresentationRate;
noiseMin = round(minNoiseDuration/vizPresentationRate); %in frames
signalMin = round(minSignalDuration/vizPresentationRate); %in frames
noisePDF = discrete_bounded_hazard_rate(expLambda, maxNoiseFrames);
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

noise1Frames = Shuffle(noiseFrames + noiseMin);
noise2Frames = Shuffle(noise1Frames);
signal1Frames = Shuffle(noiseFrames + signalMin);
if noisePeriods==0
    nFrames = noNoiseTrialDuration / vizPresentationRate;
else
    nFrames = max(signal1Frames) + 2*maxNoiseFrames + signalMin + secondSignalMin;
end
% ensure an event amount of frames for computational convenience
if mod(nFrames,2)>0
    nFrames=nFrames+1;
end
signal2Frames = nFrames - (noise1Frames + noise2Frames + signal1Frames);
trialDuration = nFrames*vizPresentationRate;

%% create arrays to hold values
choices = zeros(nTrial, 2); %raw choice; forced choice (state of accumulator)
RTs = zeros(nTrial, 1);
startPoints = zeros(nTrial, 2); %DV start point in first and second noise periods
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

% create memory evidence stream
memoryStream = (binornd(1, cue, [nFrames, nTrial])*2-1) + normrnd(0,1, [nFrames, nTrial]);
% define analytic solution for memory
% expectedMemValue = repelem([1:memoryThinning]', round(nFrames/memoryThinning))*cue;

% define viz evidence changepoints 
signal1Onsets = noise1Frames + 1;
noise2Onsets = noise1Frames + signal1Frames + 1;
signal2Onsets = noise1Frames + signal1Frames + noise2Frames + 1;
% initialize analytic solution matrix for vision
% expectedVizValue = zeros(nFrames, nTrial);

% create visual evidence stream
if noisePeriods==0
    flickerStream = (binornd(1,coherence, [nFrames,nTrial])*2-1) + normrnd(0,1,[nFrames,nTrial]);
    flickerStream(:, congruentTrials+1:nTrial) = -1*flickerStream(:, congruentTrials+1:nTrial);
else
    % make noise matrices
    if strcmp(flickerPaddingValue,'zero')==1
        flickerNoise = zeros(nFrames, nTrial);
    elseif strcmp(flickerPaddingValue, 'gaussian')==1
        flickerNoise = normrnd(0,1, [nFrames, nTrial]);
    end

    if flickerAdditiveNoise==1 && strcmp(flickerAdditiveNoiseValue, 'gaussian')==1
        flickerSampleNoise = normrnd(0,1, [nFrames,nTrial]);
    end

    % initialize stream
    if flickerNoisePadding==0
        flickerStream = ones(nFrames/2, nTrial);
    else
        flickerStream = repmat([1; NaN], nFrames/2, nTrial);
        noiseBool = isnan(flickerStream);
        flickerStream(noiseBool) = flickerNoise(noiseBool);
    end

    % populate according to coherence & noise frames
    for trial=1:nTrial
        flickerStream(1:noise1Frames(trial), trial)= flickerNoise(1:noise1Frames(trial), trial);
        flickerStream(noise2Onsets(trial):signal1Onsets(trial)-1, trial) = flickerNoise(noise2Onsets(trial):signal1Onsets(trial)-1, trial);
        imgIdx = find(flickerStream(:, trial)==1);
        nImgFrames = length(imgIdx);
        nTargetFrames = ceil(nImgFrames*coherence);
        targetIdx = randsample(imgIdx, nTargetFrames);
        lureIdx = setdiff(imgIdx, targetIdx);
        if trial <= congruentTrials
            target=1;
        else
            target=-1;
        end
        if flickerAdditiveNoise==0
            flickerStream(targetIdx, trial) = target;
            flickerStream(lureIdx, trial) = -target;
        else
            flickerStream(targetIdx, trial) = target + flickerSampleNoise(targetIdx, trial);
            flickerStream(lureIdx, trial) = -target + flickerSampleNoise(lureIdx, trial);
        end

        % compute analytic solution for each trial
        %sig1EV = (repelem([1:ceil(signal1Frames(trial)/2)]', 2)*coherence);
        %sig2EV = (repelem([1:ceil(signal2Frames(trial)/2)]', 2)*coherence) + sig1EV(signal1Frames(trial));
        %expectedVizValue(noise1Frames(trial)+1:noise2Onsets(trial)-1, trial) = sig1EV(1:signal1Frames(trial));
        %expectedVizValue(noise2Onsets(trial):signal2Onsets(trial)-1, trial) = sig1EV(signal1Frames(trial));
        %expectedVizValue(signal2Onsets(trial):nFrames, trial) = sig2EV(1:signal2Frames(trial));
    end
end


%% run simulation
tic
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
            visionPDF = betapdf(p, alphaVis, betaVis) ./ sum(betapdf(p, alphaVis, betaVis));
            % compute precision as inverse belief-weighted entropy
            visionPrecisions(frame, trial) = 1/sum(visionPDF .* entropy');
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
            memoryPrecisions(frame, trial) = 1/sum(memoryPDF .* entropy');
        else
            memoryPrecisions(frame,trial) = memoryPrecisions(frame-1,trial);
        end

        % store counter values
        counters(frame, alphaMem_idx, trial) = alphaMem;
        counters(frame, betaMem_idx, trial) = betaMem;

        % compute relative precision evidence weights
        visionDriftRate = visionPrecisions(frame,trial) / (visionPrecisions(frame,trial) + memoryPrecisions(frame,trial));
        memoryDriftRate = memoryPrecisions(frame,trial) / (visionPrecisions(frame,trial) + memoryPrecisions(frame,trial));
        visionDrift(frame, trial) = visionDriftRate;
        memoryDrift(frame, trial) = memoryDriftRate;

        % compute time-varying relative precision-weighted decision variable
        memorySample = memoryStream(frame, trial);
        visualSample = flickerStream(frame, trial);

        if frame == 1
            visionAccumulator(frame, trial) = visualSample * visionDriftRate;
            if mod(frame, memoryThinning)==1
                memoryAccumulator(frame, trial) = memorySample * memoryDriftRate;
                decisionVariable(frame, trial) = memorySample*memoryDriftRate + visualSample * visionDriftRate;
            else
                decisionVariable(frame, trial) = visualSample*visionDriftRate;
            end

        else % for frames > 1
            visionAccumulator(frame, trial) = visionAccumulator(frame-1, trial) + visualSample * visionDriftRate;
            if mod(frame, memoryThinning) == 1
                memoryAccumulator(frame, trial) = memoryAccumulator(frame-1, trial) + memorySample * memoryDriftRate;
                decisionVariable(frame, trial) = decisionVariable(frame-1, trial) + memorySample*memoryDriftRate + visualSample * visionDriftRate;
            else
                memoryAccumulator(frame, trial) = memoryAccumulator(frame-1, trial);
                decisionVariable(frame, trial) = decisionVariable(frame-1, trial) + visualSample*visionDriftRate;
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

    % populate startPoints matrix
    startPoints(trial, 1) = decisionVariable(noise1Frames(trial), trial);
    startPoints(trial, 2) = decisionVariable(signal2Onsets(trial)-1, trial);
end
toc

%% store results
% simulation settings
data.nTrial = nTrial;
data.subID = subID;
data.nSub = nSub;
data.trialDuration = trialDuration;
data.cue = cue;
data.coherence = coherence;
data.threshold = threshold;
data.congruent = logical(congruent);
data.memoryThinning = memoryThinning;
data.visionThinning = visionThinning;
data.nFrames = nFrames;
data.trialDuration = trialDuration;
data.vizPresentationRate = vizPresentationRate;
data.maxNoiseDuration = maxNoiseDuration;
data.minNoiseDuration = minNoiseDuration;
data.minSignalDuration = minSignalDuration;
data.noise1Frames = noise1Frames;
data.noise2Frames = noise2Frames;
data.noise2Onset = noise2Onsets;
data.signal2Onsets = signal2Onsets;
data.signal1Frames = signal1Frames;
data.signal2Frames = signal2Frames;
data.secondSignalMin = secondSignalMin;
data.flickerAdditiveNoise = flickerAdditiveNoise;
data.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
data.flickerNoisePadding = flickerNoisePadding;
data.flickerPaddingValue = flickerPaddingValue;
data.counters = counters;
data.startPoints = startPoints;

% behavior
data.choices = choices;
data.RT = RTs;

% optional frame-by-frame info
if saveEvidence==1
    data.memoryEvidence = memoryStream;
    data.visionEvidence = flickerStream;
end

if saveFlickerNoise==1
    data.flickerNoise = flickerNoise;
end

if saveAccumulators==1
    data.memoryAccumulator = memoryAccumulator;
    data.visionAccumulator = visionAccumulator;
end

if saveDV==1
    data.decisionVariable = decisionVariable;
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

%% save results

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

outfile = [outDir '/' num2str(cue) 'cue_' num2str(coherence) 'coh_' num2str(threshold) 'thresh_' num2str(memoryThinning) 'thin_' num2str(nTrial) 'trial_' 'sub' num2str(subID) '.mat'];
save(outfile, 'data')
end




