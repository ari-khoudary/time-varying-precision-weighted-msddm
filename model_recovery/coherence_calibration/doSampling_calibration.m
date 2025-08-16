function data = doSampling_calibration(config)

%% unpack params
nTrial = config.nTrial; 
subID = config.subID;
coherence = config.coherence;
threshold = config.threshold;
vizPresentationRate = config.vizPresentationRate;
noNoiseTrialDuration = config.noNoiseTrialDuration;
flickerPadding = config.flickerNoisePadding;

% reset random seed on every iteration for maximum stochasticity
rng('shuffle');
%% setup visual evidence & DVs on each trial

% number of frames on each trial
nFrames = noNoiseTrialDuration / vizPresentationRate;

if flickerPadding == 1
    % rows = timesteps, columns = trials
    visionEvidence = (binornd(1,coherence, [nFrames*2,nTrial])*2-1) + normrnd(0,1,[nFrames*2,nTrial]);
    noiseFrames = 1:2:size(visionEvidence, 1); % returns the index of every other row
    visionEvidence(noiseFrames, :) = 0;
else
    visionEvidence = (binornd(1,coherence, [nFrames*2,nTrial])*2-1) + normrnd(0,1,[nFrames*2,nTrial]);
end

% computes decision variable on each trial
decisionVariable = cumsum(visionEvidence, 1);

%% get RTs and accuracy for each level of threshold

% identify trials where DV hit either threshold
mask = abs(decisionVariable) > threshold;

% get index of first threshold passing
[~, rawChoiceRT] = max(mask, [], 1); 
rawChoiceRT = rawChoiceRT'; % transform from row to column vector

% indicate which trials never made it to threshold
noChoice = ~any(mask, 1)';

% set the rawChoiceRTs for those trials as NaN
rawChoiceRT(noChoice) = NaN;

% determine if DV value at threshold crossing corresponds to correct answer (positive value)
rawChoiceAccuracy = NaN(nTrial, 1);
forcedChoiceAccuracy = NaN(nTrial, 1);

for i = 1:nTrial
    if isnan(rawChoiceRT(i))
        forcedChoiceAccuracy(i) = decisionVariable(end, i) > 0;
        rawChoiceAccuracy(i) = 0;
    else
        rawChoiceAccuracy(i) = decisionVariable(rawChoiceRT(i), i) > 0;
        forcedChoiceAccuracy(i) = rawChoiceAccuracy(i);
    end
end

%% write out to table
subID = repelem(subID, nTrial)';
coherence = repelem(coherence, nTrial)';
threshold = repelem(threshold, nTrial)';
trial = 1:nTrial;
trial = trial';

data = table(subID, trial, coherence, threshold, rawChoiceRT, rawChoiceAccuracy, noChoice, forcedChoiceAccuracy);
end

