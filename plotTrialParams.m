% plot precisions, drift rates, accumulators for single trials

trial = 1;
% unwrap structure if loading in data
% struct2vars(data);
% make dummy mem evidence variable for unweighted accumulator plot
memEvidence = zeros(nFrames, 1);
memEvidence(mod(1:nFrames, memoryThinning)==1) = ...
    memoryEvidence(mod(1:nFrames, memoryThinning)==1, trial);

% set up figure
time = 1:nFrames;
alpha = 0.15;
t = tiledlayout('flow');

% estimated precision
nexttile;
hold on
plot(time, memoryPrecisions(:,trial));
plot(time, visionPrecisions(:,trial));
xregion(0, noise1Frames(trial), FaceAlpha=alpha);
xregion(noise2Onsets(trial), signal2Onsets(trial), FaceAlpha=alpha);
legend({'memory precision', 'vision precision'});
title('estimated precision');

% drift rates
nexttile;
hold on
plot(time, memoryDrifts(:, trial));
plot(time, visionDrifts(:, trial));
plot(time, visionDrifts(:,trial) + memoryDrifts(:,trial));
xregion(0, noise1Frames(trial), FaceAlpha=alpha);
xregion(noise2Onsets(trial), signal2Onsets(trial), FaceAlpha=alpha);
legend({'memory drift', 'vision drift', 'DV drift (mem+viz)'});
title('drift rates')

% accumulators
nexttile;
hold on
plot(time, memoryAccumulator(:,trial));
plot(time, visionAccumulator(:, trial));
plot(time, decisionVariable(:, trial));
xregion(0, noise1Frames(trial), FaceAlpha=alpha);
xregion(noise2Onsets(trial), signal2Onsets(trial), FaceAlpha=alpha);
legend({'memory accumulator', 'vision accumulator', 'DV'});
title('reliability-weighted accumulators')

% cumsum
nexttile;
hold on
plot(time, cumsum(memEvidence));
plot(time, cumsum(visionEvidence(:, trial)));
plot(time, decisionVariable(:, trial));
xregion(0, noise1Frames(trial), FaceAlpha=alpha);
xregion(noise2Onsets(trial), signal2Onsets(trial), FaceAlpha=alpha);
legend({'memory', 'vision', 'DV (reliability-weighted)'});
title('unweighted cumulative sum')

title(t, ['example trial (#' num2str(trial) '): ' num2str(cue) ' cue, ' num2str(coherence) ' coh, 1:' num2str(memoryThinning) ' mem:viz'])
subtitle(t, ['gray bars = noise periods; flickerAdditiveNoise = ' num2str(flickerAdditiveNoise)]);
