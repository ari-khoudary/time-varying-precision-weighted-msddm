function [choices, memoryEvidence, visualEvidence, fullEvidence, ... 
    visualDriftRate, memoryDriftRate] = doSampling(nSub, nTrial, cueLevel, coherenceLevel, congruent)

% define sampling windows
nSampMemory = (750+500)/50; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define sampling bounds
threshold = 3;
memoryStartingPoint  = 0;  
visualStartingPoint = 0;

% create variables to store values
memoryEvidence = zeros(nSub, nTrial, nSampMemory+nSampVisual); 
visualEvidence = zeros(nSub, nTrial, nSampVisual);
fullEvidence = zeros(nSub, nTrial, nSampVisual);
choices = zeros(nSub,nTrial);
RT = zeros(nSub, nTrial);
memoryDriftRates = visualEvidence;
visualDriftRates = visualEvidence;

% initialize other variables
memoryPrecision = 1/computeEntropy(cueLevel);
framesTargetA = 0;
framesTargetB = 0;

% plot "trial traces"
figure;
hold on;

for subj=1:nSub
    for trial=1:nTrial
        for t=1:nSampVisual
            for i=1:nSampMemory
                % generate & accumulate memory samples
                 memoryDriftRate = memoryPrecision / (memoryPrecision + (1/computeEntropy(0.65)));
                 memoryEvidence(subj, trial, i) = randn(1)*memoryStartingPoint + ... 
                        cumsum(memoryDriftRate*(binornd(1,cueLevel)*2-1) + normrnd(0,1)); 
            end

            % generate memory & visual samples during flicker stream
            if congruent
                memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                visualSample = (binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            else
                memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                visualSample = -(binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            end

             % add visual sample to visual evidence
             if t==1
                 visualEvidence(subj, trial, t) = visualSample;
             else
                  visualEvidence(subj, trial, t) = visualEvidence(subj, trial, t-1) + visualSample;
             end

             % add memory sample to memory evidnece
             if t==1
                memoryEvidence(subj, trial, nSampMemory+t) = memoryEvidence(subj, trial, nSampMemory) + memorySample;
             else
                 memoryEvidence(subj, trial, nSampMemory+t) = memoryEvidence(subj, trial, t-1) + memorySample;
             end

             % update target counters
            if visualSample > 0 
                framesTargetA = 1 + framesTargetA;
            else
                framesTargetB = 1 + framesTargetB;
            end

            % update target probability
             probTargetA = framesTargetA / (framesTargetA + framesTargetB);

             % update visual precision
             visualPrecision = 1/computeEntropy(probTargetA);

             % compute drift rates as relative evidence precisions
             memoryDriftRate = memoryPrecision / (visualPrecision + memoryPrecision);
             visualDriftRate = visualPrecision / (visualPrecision + memoryPrecision);

             % store values
             memoryDriftRates(subj, trial, t) = memoryDriftRate;
             visualDriftRates(subj, trial, t) = visualDriftRate;

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
        RT(subj, trial) = boundaryIdx;

        % populate choice matrix accordingly
        if fullEvidence(subj,trial,boundaryIdx) > threshold
            choices(subj,trial) = 1;
        else
            choices(subj,trial) = 0;
        end

          %plot
          plot(squeeze(fullEvidence(subj,trial, :)));
         drawnow;
    end
end

outfile = sprintf('%.1fcue_%.1fcoh_%icong_%isubs_%itrials.mat', cueLevel, coherenceLevel, congruent, nSub, nTrial);
save(outfile)
end










