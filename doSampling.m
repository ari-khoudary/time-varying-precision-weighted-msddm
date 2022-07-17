function [choices, RT] = doSampling(order, cueLevel, anticipatedCoherence, coherenceLevel, congruent, save, makePlot)

% simulates an ideal observer in a cued perceptual decision task
% make csv output & plots optional arguments

%% Initialize variables
% define sampling windows
nSampMemory = (750+500)/125; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define number of subjects & trials
nSub = 1;
nTrial = 10;

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

%% Run simulation
for subj=1:nSub
    for trial=1:nTrial
        % reset counters on each trial
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
              if order==2
                alphaMem = alphaMem + memSampsTargetA;
                betaMem = betaMem + memSampsTargetB;
                memoryPrecision = 1/betaVar(alphaMem, betaMem);
                memoryDriftRate = memoryPrecision / (memoryPrecision + (1/computeEntropy(anticipatedCoherence)));
              
              else % first-order precision
                memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);
                if i == 1
                  memoryPrecision = 1/computeEntropy(cueLevel);
                else
                  memoryPrecision = 1/computeEntropy(memProbTargetA);
                end
              memoryDriftRate = memoryPrecision / (memoryPrecision + (1/computeEntropy(anticipatedCoherence)));
              end

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
                if ~mod(t, 4) % only draw a new memory sample every 4 timepoints
                    memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                end
                visualSample = (binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            else
                if ~mod(t,4)
                    memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                end
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
            
            if ~mod(t,4)
                if memorySample > 0 
                    memSampsTargetA = 1 + memSampsTargetA;
                else
                    memSampsTargetB = 1 + memSampsTargetB;
                end
            end

            % update memory probability
            if order ==  2
                if ~mod(t,4)
                    alphaMem = alphaMem + memSampsTargetA;
                    betaMem = betaMem + memSampsTargetB;
                    memoryPrecision = 1/betaVar(alphaMem, betaMem);
                end
            % and precision-weighted drift rate
                alphaVis = alphaVis + framesTargetA;
                betaVis = betaVis + framesTargetB;
                visualPrecision = 1/betaVar(alphaVis, betaVis);

            else % first-order analogs
                 % memory probability & precision
                 if ~mod(t,4)
                    memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);
                    memoryPrecision = 1/computeEntropy(memProbTargetA);
                 end
                 
                 % vision probability & precision
                 probTargetA = framesTargetA / (framesTargetA + framesTargetB);
                 if t==1
                 visualPrecision = 1/computeEntropy(0.65);
                else
                visualPrecision = 1/computeEntropy(probTargetA);
                 end
            end

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
        if congruent 
            if fullEvidence(subj,trial,boundaryIdx) > threshold
                choices(subj,trial) = 1;
            else
                choices(subj,trial) = 0;
            end
        else
            if fullEvidence(subj,trial,boundaryIdx) < -threshold
                choices(subj,trial) = 1;
            else
                choices(subj,trial) = 0;
            end
        end
    end
end

% save results
if save
    if order == 2
        outfile = sprintf('bayes_results/%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    else
        outfile = sprintf('entropy_results/%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    end
    save(outfile)

% write csv with choices & RTs so i can plot those in R
    accuracy = choices';
    RT = RT';
    cueLevel = repmat(cueLevel, [nTrial 1]);
    anticipatedCoherence = repmat(anticipatedCoherence, [nTrial 1]);
    coherenceLevel = repmat(coherenceLevel, [nTrial 1]);
    congruent = repmat(congruent, [nTrial 1]);
    behav_csv = table(accuracy, RT, cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    if order == 2 
        filename = ['bayes_results/behav_csv/' char(extractBetween(outfile, '/', '.mat')) '.csv'];
    else
        filename = ['entropy_results/behav_csv/' char(extractBetween(outfile, '/', '.mat')) '.csv'];
    end
    writetable(behav_csv, filename);
end

%% Plot results
% trial traces
if makePlot
    addX = NaN(nSampMemory,1);
    fig=figure; 
    for a=1:10
        subplot(2, 5, a)
        hold on;
        plot(squeeze(memoryEvidence(1,a,:)), 'LineWidth',1.5)
        plot([addX;squeeze(visualEvidence(1, a, :))], 'LineWidth',1.5)
        plot([addX;squeeze(fullEvidence(1,a,:))], 'LineWidth',1.5)
        plot([1,140],[0,0],'k')
        xline(nSampMemory, 'k--')
        yline([threshold -threshold])
        string2 = sprintf('trial %i', a);
        title(string2);
    end
    if congruent
        string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, congruent', cueLevel, anticipatedCoherence, coherenceLevel);
    else
        string = sprintf('%.2f cue, %.2f anticipated coherence, %.2f actual coherence, incongruent', cueLevel, anticipatedCoherence, coherenceLevel);
    end
    if order==1
        sgtitle(sprintf(['first-order model\n' string]));
    else
        sgtitle(sprintf(['second-order model\n' string]));
    end
    h=legend({'memory','visual','combined'},'FontSize',8, 'Orientation', 'horizontal');
    set(h, 'Position', [0.65 0.46 0.25 0.025]);
    plots=axes(fig, 'visible', 'off');
    plots.XLabel.Visible='on';
    plots.YLabel.Visible='on';
    plots.Title.Visible='on';
    xlabel(plots, 'time (a.u.)');
    ylabel(plots, 'evidence (a.u.)');
    if save
        outfig = char(regexp(outfile, 'b.*cong', 'match'));
        saveas(gcf, [outfig '_traces.png']);
    end
     
    % drifts
    fig=figure; 
    for b=1:10
        subplot(2, 5, b)
        hold on;
        plot(squeeze(memoryDriftRates(1,b,:)), 'LineWidth',1.5)
        plot([addX;squeeze(visualDriftRates(1, b, :))], 'LineWidth',1.5)
        plot([1,140],[0,0],'k')
        xline(nSampMemory, 'k--')
        string2 = sprintf('trial %i', b);
        title(string2);
    end
    h=legend({'memory','visual'},'Orientation','horizontal');
    set(h, 'Position', [0.65 0.46 0.25 0.025]);
    if order==1
        sgtitle(sprintf(['first-order model\n' string]));
    else
        sgtitle(sprintf(['second-order model\n' string]));
    end
    plots=axes(fig, 'visible', 'off');
    plots.XLabel.Visible='on';
    plots.YLabel.Visible='on';
    plots.Title.Visible='on';
    xlabel(plots, 'time (a.u.)');
    ylabel(plots, 'drift rate (a.u.)');
    if save
        saveas(gcf, [outfig '_drifts.png']);
    end
end
end










