function data = doSampling_calibration(config)

%% unpack params
nTrial = config.nTrial; 
subID = config.subID;
coherence = config.coherence;
thresholds = config.threshold;
vizPresentationRate = config.vizPresentationRate;
noNoiseTrialDuration = config.noNoiseTrialDuration;
flickerPadding = config.flickerNoisePadding;

% where do you want to save the results?
outDir = config.outDir;

%% setup
nFrames = noNoiseTrialDuration / vizPresentationRate;

if flickerPadding == 1
    visionEvidence = (binornd(1,coherence, [nFrames*2,nTrial])*2-1) + normrnd(0,1,[nFrames*2,nTrial]);
    for trial = 1:width(visionEvidence)
        for t = 1:height(visionEvidence)
            if mod(t, 2) ~= 0
                visionEvidence(t, trial) = 0;
            end
        end
    end
else
    visionEvidence = (binornd(1,coherence, [nFrames*2,nTrial])*2-1) + normrnd(0,1,[nFrames*2,nTrial]);
end

decisionVariable = cumsum(visionEvidence, 1);

data = table('Size', [nTrial length(thresholds*2)]);

%% get RTs and accuracy

for i = 1:length(thresholds)
    
    threshold = thresholds(i);
    mask = decisionVariable > threshold;

    [~, rawChoiceRT] = max(mask, [], 1);
    rawChoiceAccuracy = any(mask, 1);   
    forcedChoiceRT(~rawChoiceAccuracy) = height(decisionVariable);
    rawChoiceRT(~rawChoiceAccuracy) = NaN;
    forcedChoiceAccuracy = rawChoiceAccuracy | (decisionVariable(end, :) > threshold);

end

