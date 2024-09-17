function data = doSampling_noMemory(config)

%% unpack params
nTrial = config.nTrial; %per cue
nSub = config.nSub;
subID = config.subID;
coherence = config.coherence;
threshold = config.threshold;
visionThinning = config.visionThinning;
vizPresentationRate = config.vizPresentationRate;
noisePeriods = config.noisePeriods;
expLambda = config.expLambda;
maxNoiseDuration = config.maxNoiseDuration;
minNoiseDuration = config.minNoiseDuration;
minSignalDuration = config.minSignalDuration;
secondSignalMin = config.secondSignalMin;
noNoiseTrialDuration = config.noNoiseTrialDuration;
flickerAdditiveNoise = config.flickerAdditiveNoise;
flickerAdditiveNoiseValue = config.flickerAdditiveNoiseValue;
flickerPadding = config.flickerNoisePadding;
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
% make noise durations
if noisePeriods == 1
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

    noise1Frames = noiseFrames + noiseMin;
    noise1Frames = noise1Frames(randperm(length(noise1Frames)));
    noise2Frames = noise1Frames(randperm(length(noise1Frames)));
    signal1Frames = noiseFrames + signalMin;
    signal1Frames = signal1Frames(randperm(length(signal1Frames)));
end

if noisePeriods==0
    nFrames = noNoiseTrialDuration / vizPresentationRate;
else
    nFrames = round(max(signal1Frames) + 2*maxNoiseFrames + signalMin + secondSignalMin);
end
% ensure an event amount of frames for computational convenience
if mod(nFrames,2)>0
    nFrames=nFrames+1;
end
if noisePeriods==1
    signal2Frames = nFrames - (noise1Frames + noise2Frames + signal1Frames);
end
trialDuration = nFrames*vizPresentationRate;

%% create arrays to hold values
choices = zeros(nTrial, 2); %raw choice; forced choice (state of accumulator)
RTs = zeros(nTrial, 1);
startPoints = zeros(nTrial, 2); %DV start point in first and second noise periods
accumulatorPosition = zeros(nTrial, 1); % value of vis accumulator at RT
visionAccumulator = zeros(nFrames, nTrial);
decisionVariable = zeros(nFrames, nTrial);
visionPrecisions = zeros(nFrames, nTrial);
visionDrift = zeros(nFrames, nTrial);
counters = zeros(nFrames, 4, nTrial);

% create index variables
alphaVis_idx = 3;
betaVis_idx = 4;

%% set up
% x vector for PDFs
p = 0.01:.01:0.99;

% create visual evidence stream
if noisePeriods==0
    visionEvidence = (binornd(1,coherence, [nFrames,nTrial])*2-1) + normrnd(0,1,[nFrames,nTrial]);
    for i = 1:size(visionEvidence, 2)
        for j = 1:size(visionEvidence, 1)
            if mod(j, 2) == 0
                visionEvidence(j, i) = 0; % incorporate inter-stimulus masks
            end
        end
    end
else
    % define viz evidence changepoints
    signal1Onsets = noise1Frames + 1;
    noise2Onsets = noise1Frames + signal1Frames + 1;
    signal2Onsets = noise1Frames + signal1Frames + noise2Frames + 1;
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
    if flickerPadding==0
        visionEvidence = ones(nFrames/2, nTrial);
    else
        visionEvidence = repmat([1; NaN], nFrames/2, nTrial);
        noiseBool = isnan(visionEvidence);
        visionEvidence(noiseBool) = flickerNoise(noiseBool);
    end

    % populate according to coherence & noise frames
    for trial=1:nTrial
        visionEvidence(1:noise1Frames(trial), trial)= flickerNoise(1:noise1Frames(trial), trial);
        visionEvidence(noise2Onsets(trial):signal2Onsets(trial)-1, trial) = flickerNoise(noise2Onsets(trial):signal2Onsets(trial)-1, trial);
        imgIdx = find(visionEvidence(:, trial)==1);
        nImgFrames = length(imgIdx);
        nTargetFrames = ceil(nImgFrames*coherence);
        targetIdx = randsample(imgIdx, nTargetFrames);
        lureIdx = setdiff(imgIdx, targetIdx);
        target = 1;
        if flickerAdditiveNoise==0
            visionEvidence(targetIdx, trial) = target;
            visionEvidence(lureIdx, trial) = -target;
        else
            visionEvidence(targetIdx, trial) = target + flickerSampleNoise(targetIdx, trial);
            visionEvidence(lureIdx, trial) = -target + flickerSampleNoise(lureIdx, trial);
        end
    end
end


%% run simulation
for trial=1:nTrial
    % initialize counters
    alphaVis = 1;
    betaVis = 1;

    % compute time-varying drift rate
    for frame=1:nFrames

        if mod(frame, visionThinning)==0
            % update visual counters
            if visionEvidence(frame, trial) > 0
                alphaVis = alphaVis + 1;
            elseif visionEvidence(frame, trial) < 0
                betaVis = betaVis + 1;
            end

            % compute normalized beta pdf
            visionPDF = betapdf(p, alphaVis, betaVis) ./ sum(betapdf(p, alphaVis, betaVis));
            % compute precision
            visionEntropy = computeEntropy(visionPDF);
            visionPrecisions(frame, trial) = 1/visionEntropy;
        end

        % store counter values
        counters(frame, alphaVis_idx, trial) = alphaVis;
        counters(frame, betaVis_idx, trial) = betaVis;

        % compute drift rate as visual evidence reliability
        visionDriftRate = visionPrecisions(frame, trial);

        % compute time-varying relative precision-weighted decision variable
        visualSample = visionEvidence(frame, trial);

        if frame == 1
            visionAccumulator(frame, trial) = visualSample * visionDriftRate;
            % decisionVariable(frame, trial) = bias + visualSample*visionDriftRate;
        else % for frames > 1
            visionAccumulator(frame, trial) = visionAccumulator(frame-1, trial) + visualSample * visionDriftRate;
            decisionVariable(frame, trial) = decisionVariable(frame-1, trial) + visualSample*visionDriftRate;
        end

        % make decision
        % find point at which evidence crosses threshold -- USING JUST VIZ ACCUMULATOR 
        boundaryIdx = find(abs(visionAccumulator(:, trial))>threshold, 1);
        if isempty(boundaryIdx)
            boundaryIdx = nFrames;
        end

        % store crossing point as RT
        RTs(trial) = boundaryIdx;

        % populate choice matrix accordingly - first column is "raw"
        % responses, second is "forced"
        if visionAccumulator(boundaryIdx, trial) > threshold
            choices(trial, 1) = 1;
            choices(trial, 2) = 1;
        elseif boundaryIdx==nFrames
            choices(trial, 1) = NaN;
            if visionAccumulator(boundaryIdx, trial) > 0
                choices(trial, 2) = 1;
            else
                choices(trial, 2) = 0;
            end
        else % if wrong bound is reached
            choices(trial, 1) = 0;
            choices(trial, 2) = 0;
        end

        % populate startPoints matrix
        if noisePeriods==1
            startPoints(trial, 1) = decisionVariable(noise1Frames(trial), trial);
            startPoints(trial, 2) = decisionVariable(signal2Onsets(trial)-1, trial);
        else
            startPoints(trial, 1) = NaN;
            startPoints(trial, 2) = NaN;
        end

        % store terminal values of non-DV accumulators
        accumulatorPosition(trial) = visionAccumulator(boundaryIdx, trial);
    end
end

    %% store results
    % simulation settings
    data.nTrial = nTrial;
    data.subID = subID;
    data.nSub = nSub;
    data.trialDuration = trialDuration;
    data.coherence = coherence;
    data.threshold = threshold;
    data.accumulatorPosition = accumulatorPosition;
    data.visionThinning = visionThinning;
    data.nFrames = nFrames;
    data.trialDuration = trialDuration;
    data.vizPresentationRate = vizPresentationRate;
    data.maxNoiseDuration = maxNoiseDuration;
    data.minNoiseDuration = minNoiseDuration;
    data.minSignalDuration = minSignalDuration;
    data.noisePeriods = noisePeriods;
    if noisePeriods==1
        data.noise1Frames = noise1Frames;
        data.noise2Frames = noise2Frames;
        data.signal1Onsets = signal1Onsets;
        data.noise2Onsets = noise2Onsets;
        data.signal2Onsets = signal2Onsets;
        data.signal1Frames = signal1Frames;
        data.signal2Frames = signal2Frames;
        data.secondSignalMin = secondSignalMin;
    end
    data.flickerAdditiveNoise = flickerAdditiveNoise;
    data.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
    data.flickerNoisePadding = flickerPadding;
    data.flickerPaddingValue = flickerPaddingValue;
    data.counters = counters;
    data.startPoints = startPoints;

    % behavior
    data.choices = choices;
    data.RT = RTs;

    % optional frame-by-frame info
    if saveEvidence==1
        data.visionEvidence = visionEvidence;
    end

    if saveFlickerNoise==1
        data.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
    end

    if saveAccumulators==1
        data.visionAccumulator = visionAccumulator;
    end

    if saveDV==1
        data.decisionVariable = decisionVariable;
    end

    if saveCounters==1
        data.counters=counters;
    end

    if savePrecisions==1
        data.visionPrecisions = visionPrecisions;
    end

    if saveDrifts==1
        data.visionDrifts = visionDrift;
    end

    %% save results

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    outfile = [outDir 'noMemory_' num2str(coherence) 'coh_' num2str(threshold) 'thresh_' num2str(nTrial) 'trial_' 'sub' num2str(subID) '.mat'];
    save(outfile, 'data')
end



