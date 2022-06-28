function [choices, RT] = doSampling(cueLevel, coherenceLevel, congruent, anticipatedCoherence)

% define sampling windows
nSampMemory = (750+500)/50; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define number of subjects & trials
nSub = 1;
nTrial = 10;

% define sampling bounds
threshold = 10;
memoryStartingPoint = 0;  
visualStartingPoint = 0;

% create variables to store values
memoryEvidence = zeros(nSub, nTrial, nSampMemory+nSampVisual); 
visualEvidence = zeros(nSub, nTrial, nSampVisual);
fullEvidence = zeros(nSub, nTrial, nSampVisual);
choices = zeros(nSub,nTrial);
RT = zeros(nSub, nTrial);
memoryDriftRates = memoryEvidence;
visualDriftRates = visualEvidence;
memoryPrecisions = memoryEvidence;
visualPrecisions = visualEvidence;

% initialize other variables
betaX = 0:.01:1;

% run simulation
for subj=1:nSub
    for trial=1:nTrial
        % reset all counters at the start of each trial
        alphaMem = 1;
        betaMem = 1;
        alphaVis = 1;
        betaVis = 1;
        framesTargetA = 0;
        framesTargetB = 0;
        memSampsTargetA = 0;
        memSampsTargetB = 0;

        % start with memory retrieval period
        for i=1:nSampMemory
            % generate memory sample
              memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
              % store outcome of sample
              if memorySample> 0
                 memSampsTargetA = 1 + memSampsTargetA;
              else
                  memSampsTargetB = 1 + memSampsTargetB;
              end

              % compute precision-weighted drift rate
              alphaMem = alphaMem + memSampsTargetA;
              betaMem = betaMem + memSampsTargetB;
              %memProbTargetA = betapdf(betaX, alphaMem, betaMem);
              memoryRetrievalPrecision = 1/betaVar(alphaMem, betaMem);

              if exist('anticipatedCoherence', 'var')
                memoryDriftRate = memoryRetrievalPrecision / (memoryRetrievalPrecision + (1/computeEntropy(anticipatedCoherence)));
              else
                memoryDriftRate = memoryRetrievalPrecision / (memoryRetrievalPrecision + 1/var(betapdf(betaX, alphaVis, betaVis)));
              end

              % store values
              memoryPrecisions(subj, trial, i) = memoryRetrievalPrecision;
              memoryDriftRates(subj, trial, i) = memoryDriftRate;

              % accumulate evidence
              if i == 1
                 memoryEvidence(subj, trial, i) = randn(1)+memoryStartingPoint + memoryDriftRate*memorySample;
             else
                memoryEvidence(subj, trial, i) = memoryEvidence(subj, trial, i-1) + memoryDriftRate*memorySample;
              end
        end

        % then begin parallel visual & memory sampling
        for t=1:nSampVisual

            % generate memory & visual samples
            if congruent
                memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                visualSample = (binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            else
                memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                visualSample = -(binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            end

             % add samples to evidence matrices
             if t==1
                 visualEvidence(subj, trial, t) = visualSample;
                 memoryEvidence(subj, trial, nSampMemory+t) = memoryEvidence(subj, trial, nSampMemory) + memorySample;
             else
                  visualEvidence(subj, trial, t) = visualEvidence(subj, trial, t-1) + visualSample;
                  memoryEvidence(subj, trial, nSampMemory+t) = memoryEvidence(subj, trial, nSampMemory+(t-1)) + memorySample;
             end

             % update target counters
            if visualSample > 0 
                framesTargetA = 1 + framesTargetA;
            else
                framesTargetB = 1 + framesTargetB;
            end

            if memorySample > 0 
                memSampsTargetA = 1 + memSampsTargetA;
            else
                memSampsTargetB = 1 + memSampsTargetB;
            end

            % update memory probability
            alphaMem = alphaMem + memSampsTargetA;
            betaMem = betaMem + memSampsTargetB;
            %memProbTargetA = betapdf(betaX, alphaMem, betaMem);
            memoryPrecision = 1/betaVar(alphaMem, betaMem);

            % compute precision-weighted drift rate
            alphaVis = alphaVis + framesTargetA;
            betaVis = betaVis + framesTargetB;
            %probTargetA = betapdf(betaX, alphaVis, betaVis);
            visualPrecision = 1/betaVar(alphaVis, betaVis);


             % compute drift rates as relative evidence precisions
             memoryDriftRate = memoryPrecision / (visualPrecision + memoryPrecision);
             visualDriftRate = visualPrecision / (visualPrecision + memoryPrecision);

             % store values
             memoryDriftRates(subj, trial, nSampMemory+t) = memoryDriftRate;
             visualDriftRates(subj, trial, t) = visualDriftRate;
             memoryPrecisions(subj, trial, nSampMemory+t) = memoryPrecision;
             visualPrecisions(subj, trial, t) = visualPrecision;

             % accumulate precision-weighted additive evidence
             if t==1
                fullEvidence(subj, trial, t) = randn(1)+visualStartingPoint + memorySample*memoryDriftRate + visualSample*visualDriftRate;
             else
                 fullEvidence(subj, trial, t) = cumsum(fullEvidence(subj, trial, t-1) + ...
                     memorySample*memoryDriftRate + visualSample*visualDriftRate);
             end
        end
        
        % find point at which evidence crosses threshold
        boundaryIdx = find(abs(fullEvidence(subj,trial,:))>threshold, 1);
        if isempty(boundaryIdx)
            boundaryIdx = nSampVisual;
        end
        RT(subj, trial) = boundaryIdx;

        % populate choice matrix accordingly
        if fullEvidence(subj,trial,boundaryIdx) > threshold
            choices(subj,trial) = 1;
        else
            choices(subj,trial) = 0;
        end
    end
end

if exist('anticipatedCoherence', 'var')
    outfile = sprintf('results/%.2fcue_%.2fcoh_%icong_%isubs_%itrials_%.2fantCoh_bayesUpdate.mat', cueLevel, coherenceLevel, congruent, nSub, nTrial,anticipatedCoherence);
else
    outfile = sprintf('results/%.2fcue_%.2fcoh_%icong_%isubs_%itrials_0antCoh_bayesUpdate.mat', cueLevel, coherenceLevel, congruent, nSub, nTrial);
end
save(outfile)
end










