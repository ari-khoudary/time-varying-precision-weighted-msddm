%% Initialize variables
tic
% define task settings
writeCSV = 0;
plotResults = 0;
nSub = 1;
nTrial = 100; %per cue
trialDuration = 9; %in seconds; maximum possible duration
cueLevels = [0.7];
coherenceLevels = [0.52];
thinningLevels = [12];
thresholdLevels = [9];
visualThinning = 1;
ifi = 1/60;
nFrames = trialDuration / ifi;
maxNoiseDuration = 2; %in seconds
minNoiseDuration = 1.25; %in seconds
minSignalDuration = 0.375; %in seconds
noiseFrames = zeros(nTrial, 1);
maxNoiseFrames = maxNoiseDuration/ifi;
noiseMin = round(minNoiseDuration/ifi); %in frames
signalMin = round(minSignalDuration/ifi); %in frames


% create variables to store values
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

% precompute entropy distribution
p = 0.05:.1:0.95;
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

% create index variables
alphaMem_idx = 1;
betaMem_idx = 2;
alphaVis_idx = 3;
betaVis_idx = 4;
memory_idx = 1;
visual_idx = 2;
dv_idx = 3;

%% Run simulation
for subj=1:nSub

    if writeCSV==1
        csv_name = ['pilot_tests/round4/behav_csv/sub' num2str(subj) '.csv'];
        csvFile = fopen(csv_name, 'wt+');
        fprintf(csvFile, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s', ...
            'rawChoice', 'forcedChoice', 'RT', 'cueLevel', 'coherence', 'congruent', 'thinning', 'threshold', 'noise1Frames', 'signal1Frames', 'noise2Frames', 'subID');
    end

    % loop over coherence levels
    for coh = 1:length(coherenceLevels)
        coherence = coherenceLevels(coh);

        % then over cue levels
        for cue=1:length(cueLevels)
            cueLevel = cueLevels(cue);

            % then over threshold levels
            for thresh = 1:length(thresholdLevels)
                threshold = thresholdLevels(thresh);

                % then over thinning levels
                for thin = 1:length(thinningLevels)
                    memoryThinning = thinningLevels(thin);

                    % create evidence streams for this configuration of coherence, cue, threshold, and thinning
                    memoryStream = (binornd(1, cueLevel, [nFrames,nTrial])*2-1) + normrnd(0,1, [nFrames,nTrial]);

                    % then vision
                    noise1Frames = Shuffle(noiseFrames + noiseMin);
                    noise2Frames = Shuffle(noise1Frames);
                    signal1Frames = Shuffle(noiseFrames + signalMin);
                    noise2Onset = noise1Frames + signal1Frames;
                    congruentTrials = ceil(nTrial*cueLevel);
                    congruent = [ones(congruentTrials, 1); zeros(nTrial-congruentTrials,1)];
                    flickerStream = repmat([0; 1], nFrames/2, nTrial);

                    for trial=1:nTrial
                        flickerStream(1:noise1Frames(trial), trial)= 0;
                        flickerStream(noise2Onset(trial):noise2Onset(trial)+noise2Frames(trial), trial) = 0;
                        imgIdx = find(flickerStream(:,trial));
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

                    %% begin simulating task trials
                    for trial=1:nTrial
                        % initialize counters
                        alphaMem = 1;
                        betaMem = 1;
                        alphaVis = 1;
                        betaVis = 1;

                        % compute time-varying drift rate
                        for frame=1:nFrames

                            if mod(frame, visualThinning)==0
                                % update visual counters
                                if flickerStream(frame, trial) > 0
                                    alphaVis = alphaVis + 1;
                                elseif flickerStream(frame,trial) < 0
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
                            memoryDrift(frame,trial) = memoryDriftRate;

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
                                if decisionVariable(boundaryIdx,  trial) > 0
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

                        if writeCSV==1
                            % write trial data to csv
                            fprintf(csvFile, '\n %i, %i, %i, %.2f, %.4f, %i, %i, %i, %i, %i, %i, %i', ...
                                choices(trial,1), choices(trial,2), RTs(trial), cueLevel, coherence, congruent(trial), memoryThinning, threshold, noise1Frames(trial), signal1Frames(trial), noise2Frames(trial), subj);
                        end
                    end
                end
            end
        end
    end
end
toc

