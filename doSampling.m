function [choices, memoryEvidence, visualEvidence, fullEvidence, ... 
    visualDriftRate, memoryDriftRate] = doSampling(cueLevel, coherenceLevel, congruent)

% define sampling windows
nSampMemory = (750+500)/50; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define number of subjects & trials
nSub = 100;
nTrial = 25;

% define sampling bounds
threshold = 10;
memoryStartingPoint  = 0;  
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
framesTargetA = 0;
framesTargetB = 0;
memSampsTargetA = 0;
memSampsTargetB = 0;

% plot "trial traces"
%figure;
%hold on;

for subj=1:nSub
    for trial=1:nTrial
        for i=1:nSampMemory
            % generate memory sample
              memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
              % update target probability
              if memorySample> 0
                 memSampsTargetA = 1 + memSampsTargetA;
              else
                 memSampsTargetB = 1 + memSampsTargetB;
              end
              % compute precision-weighted drift rate
              memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);
              memoryRetrievalPrecision = 1/computeEntropy(memProbTargetA);
              memoryDriftRate = memoryRetrievalPrecision / (memoryRetrievalPrecision + (1/computeEntropy(0.65)));

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
        
        for t=1:nSampVisual

            % generate memory & visual samples during flicker stream
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
                  memoryEvidence(subj, trial, nSampMemory+t) = memoryEvidence(subj, trial, t-1) + memorySample;
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
                memSampsTargetA = 1 + memSampsTargetB;
            end

            % update target probabilities
             probTargetA = framesTargetA / (framesTargetA + framesTargetB);
             memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);

             % update precisions
             visualPrecision = 1/computeEntropy(probTargetA);
             memoryPrecision = 1/computeEntropy(memProbTargetA);

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
                fullEvidence(subj, trial, t) = visualStartingPoint + memorySample*memoryDriftRate + visualSample*visualDriftRate;
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

          %plot
          %plot(squeeze(fullEvidence(subj,trial, :)));
         %drawnow;
    end
end

outfile = sprintf('results/%.1fcue_%.1fcoh_%icong_%isubs_%itrials.mat', cueLevel, coherenceLevel, congruent, nSub, nTrial);
save(outfile)
end










