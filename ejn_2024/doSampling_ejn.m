function data = doSampling_ejn(config)

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
trialDuration = config.trialDuration;
delayPeriod = config.delayPeriod;
delayDurations = config.delayDurations;
halfNeutralTrials = config.halfNeutralTrials;

% do you want to save any frame-by-frame information?
saveEvidence = config.saveEvidence;
saveAccumulators = config.saveAccumulators;
saveDV = config.saveDV;
saveCounters = config.saveCounters;
savePrecisions = config.savePrecisions;
saveDrifts = config.saveDrifts;

% where do you want to save the results?
outDir = config.outDir;

%% translate parameters into simulation-space
rng("shuffle")
if cue==0.5 && halfNeutralTrials==1
    nTrial=nTrial/2;
end
congruentTrials = ceil(nTrial*cue);
congruent = [ones(congruentTrials, 1); zeros(nTrial-congruentTrials,1)];

% make noise durations
if delayPeriod == 1
    nFrames = round(max(delayDurations) + trialDuration) / vizPresentationRate;
else
    nFrames = trialDuration / vizPresentationRate;
end
% ensure an event amount of frames for computational convenience
if mod(nFrames,2)>0
    nFrames=nFrames+1;
end

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
memoryEvidence = (binornd(1, cue, [nFrames, nTrial])*2-1) + normrnd(0,1, [nFrames, nTrial]);

% create visual evidence stream
visionEvidence = (binornd(1,coherence, [nFrames,nTrial])*2-1) + normrnd(0,1,[nFrames,nTrial]);
visionEvidence(:, congruentTrials+1:nTrial) = -1*visionEvidence(:, congruentTrials+1:nTrial);

% compute delay periods 
if delayPeriod==1
    trialDelays = zeros(1, nTrial);
    for t = 1:nTrial
        if unifrnd(0,1) < 1/3
            trialDelays(t) = delayDurations(1)/vizPresentationRate;
        elseif unifrnd(0,1) < 2/3
            trialDelays(t) = delayDurations(2)/vizPresentationRate;
        else
            trialDelays(t) = delayDurations(3)/vizPresentationRate;
        end
        visionEvidence(1:trialDelays(t), t) = 0;
    end
end


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
            if visionEvidence(frame, trial) > 0
                alphaVis = alphaVis + 1;
            elseif visionEvidence(frame, trial) < 0
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
            if memoryEvidence(frame, trial) > 0
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
        memorySample = memoryEvidence(frame, trial);
        visualSample = visionEvidence(frame, trial);

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
    if delayDurations==1
        startPoints(trial, 1) = decisionVariable(trialDelays(trial)+1, trial);
    end
end

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
data.delayPeriod = delayPeriod;
if delayPeriod==1
    data.delayDurations = delayDurations;
    data.trialDelays = trialDelays;
else
    data.delayDurations = NaN;
    data.trialDelays = NaN;
end
data.halfNeutralTrials = halfNeutralTrials;
data.counters = counters;
data.startPoints = startPoints;

% behavior
data.choices = choices;
data.RT = RTs;

% optional frame-by-frame info
if saveEvidence==1
    data.memoryEvidence = memoryEvidence;
    data.visionEvidence = visionEvidence;
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




