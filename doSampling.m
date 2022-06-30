function [choices, RT] = doSampling(cueLevel, anticipatedCoherence, coherenceLevel, congruent)

% simulates an ideal observer in a cued perceptual decision task

% define sampling windows
nSampMemory = (750+500)/50; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define number of subjects & trials
nSub = 1;
nTrial = 5000;

% define sampling bounds
threshold = 25;
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

%%%%%%%% run simulation %%%%%%%%
for subj=1:nSub
    for trial=1:nTrial
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
            % generate & save memory sample
              memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
              if memorySample> 0
                 memSampsTargetA = 1 + memSampsTargetA;
              else
                  memSampsTargetB = 1 + memSampsTargetB;
              end

              % compute precision-weighted drift rate
              alphaMem = alphaMem + memSampsTargetA;
              betaMem = betaMem + memSampsTargetB;
              memoryPrecision = 1/betaVar(alphaMem, betaMem);
              memoryDriftRate = memoryPrecision / (memoryPrecision + (1/computeEntropy(anticipatedCoherence)));

              % store values
              memoryPrecisions(subj, trial, i) = memoryPrecision;
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
            memoryPrecision = 1/betaVar(alphaMem, betaMem);

            % compute precision-weighted drift rate
            alphaVis = alphaVis + framesTargetA;
            betaVis = betaVis + framesTargetB;
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

% save results
outfile = sprintf('bayes_results/%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent);
save(outfile)

% write csv with choices & RTs so i can plot those in R
accuracy = choices';
RT = RT';
cueLevel = repmat(cueLevel, [nTrial 1]);
anticipatedCoherence = repmat(anticipatedCoherence, [nTrial 1]);
coherenceLevel = repmat(coherenceLevel, [nTrial 1]);
congruent = repmat(congruent, [nTrial 1]);
behav_csv = table(accuracy, RT, cueLevel, anticipatedCoherence, coherenceLevel, congruent);
filename = ['bayes_results/behav_csv/' char(extractBetween(outfile, '/', '.mat')) '.csv'];
writetable(behav_csv, filename);


%%%%%%%% plot results %%%%%%%%
% trial traces
addX = NaN(nSampMemory,1);
fig=figure; 
for i=1:10
    subplot(2, 5, i)
    hold on;
    plot(squeeze(memoryEvidence(subj,i,:)), 'LineWidth',1.5)
    plot([addX;squeeze(visualEvidence(subj, i, :))], 'LineWidth',1.5)
    plot([addX;squeeze(fullEvidence(subj,i,:))], 'LineWidth',1.5)
    plot([1,140],[0,0],'k')
    xline(nSampMemory, 'k--')
    string2 = sprintf('trial %i', i);
    title(string2);
end
if congruent
    string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, congruent', cueLevel, anticipatedCoherence, coherenceLevel);
else
    string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, incongruent', cueLevel, anticipatedCoherence, coherenceLevel);
end
sgtitle(string);
h=legend({'memory','visual','combined'},'FontSize',8, 'Orientation', 'horizontal');
set(h, 'Position', [0.55 0.48 0.35 0.025]);
outfig = char(regexp(outfile, 'b.*cong', 'match'));
plots=axes(fig, 'visible', 'off');
plots.XLabel.Visible='on';
plots.YLabel.Visible='on';
plots.Title.Visible='on';
xlabel(plots, 'time (a.u.)');
ylabel(plots, 'evidence (a.u.)');
saveas(gcf, [outfig '_traces.png'])
 
% drifts
fig=figure; 
for j=1:10
    subplot(2, 5, j)
    hold on;
    plot(squeeze(memoryDriftRates(subj,j,:)), 'LineWidth',1.5)
    plot([addX;squeeze(visualDriftRates(subj, j, :))], 'LineWidth',1.5)
    plot([1,140],[0,0],'k')
    xline(nSampMemory, 'k--')
    string2 = sprintf('trial %i', i);
    title(string2);
end
h=legend({'memory','visual'},'Orientation','horizontal');
set(h, 'Position', [0.65 0.48 0.25 0.025]);
sgtitle(string);
plots=axes(fig, 'visible', 'off');
plots.XLabel.Visible='on';
plots.YLabel.Visible='on';
plots.Title.Visible='on';
xlabel(plots, 'time (a.u.)');
ylabel(plots, 'drift rate (a.u.)');
saveas(gcf, [outfig '_drifts.png'])

end










