function [choices, RT] = doSampling(order, cueLevel, anticipatedCoherence, coherenceLevel, congruent, ...
    threshold, memoryThinning, v, ...
    saveFiles, plotResults)

% multi-stage, multi-source sequential sampling model. sample weights vary
% dynamically as the relative precision of each evidence stream, updated
% after each sample.

% INPUTS
% order [1,2]: determines whether the precision will be computed in a first- or
% second-order manner. first-order precision=inverse entropy of evidence
% stream; second-order precision=inverse variance of the distribution over
% the data-generating parameter. CCN 2022 results use the second-order
% model only

% cueLevel [0:1]: sets the precision of memory evidence

% anticipatedCoherence [0:1]: sets the observer's expectation about the
% precision of upcoming visual evidence

% coherenceLevel [0:1]: sets the precision of visual evidence

% congruent [0,1]: determines whether visual evidence is congruent
% with memory evidence

% threshold [1:Inf]: bounds of accumulator

% memoryThinning [1:Inf]: how much to slow the memory sample rate relative to
% visual sampling rate during parallel sampling. memoryThinning=1 sets
% memory sample rate equal to visual sampling rate (via mod())

% v [1:Inf]: scales the size of the visual prior for the first order model,
% will be updated soon to also scale the size of the visual prior for the
% second order model. all CCN 2022 results were generated with v=0

% saveFiles [0,1]: turns off/on saving of .mat files, exporting of
% behavioral .csv files, and saving of figures

% plotResults [0,1]: turns off/on automatic plotting of evidence traces
% and drift rates for the first 10 trials of the simulation

% OUTPUTS
% choices: 1xnTrial array of optimal choices for the above-defined task
% settings

% RT: 1xnTrial array of reaction times (arbitrary units) for those choices

% MODIFIABLE SIMULATION SETTINGS
% nSampMemory [1:Inf]: number of pre-cue memory samples / duration of memory
% sampling period

% nSampVisual [1:Inf]: number of visual samples / duration of parallel sampling
% period

% nSub [1:Inf]: number of subjects. currently we do not model any individual
% differences so this parameter is irrelevant

% nTrial [1:Inf]: number of trials. if auto-plotting of results is desired, a
% minimum of 10 trials is needed

% memory or visual StartingPoint [0:Inf]: intercept of each evidence source's
% accumulator

% Code written by Ari Khoudary with help from Aaron Bornstein & Megan
% Peters

%% Initialize variables
% define sampling windows
nSampMemory = (750+500)/125; %miliseconds allotted in the experiment divided by estimate of memory sampling rate
nSampVisual = (1/30)*3000; %computer refresh rate * miliseconds allotted in the experiment

% define number of subjects & trials
nSub = 50;
nTrial = 25;

% define DDM parameters
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
memoryAlphas = memoryEvidence;
memoryBetas = memoryEvidence;
visualAlphas = visualEvidence;
visualBetas = visualEvidence;

%% Run simulation
for subj=1:nSub
        % make an array of "ghost samples" that functions as the visual prior
         ghostSamples = (binornd(1, anticipatedCoherence, [nSub, nTrial, v])*2-1) + normrnd(0,1, [nSub, nTrial, v]);
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

        % anticipatory memory sampling
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
                memoryAlphas(subj, trial, i) = alphaMem;
                memoryBetas(subj, trial, i) = betaMem;
                memoryPrecision = 1/betaVar(alphaMem, betaMem);
              
              else % first-order precision
                memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);
                if i == 1
                  memoryPrecision = 1/computeEntropy(cueLevel);
                else
                  memoryPrecision = 1/computeEntropy(memProbTargetA);
                end
              end
              
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
                if ~mod(t, memoryThinning) % only draw a new memory sample every timepoint set by memoryThinning
                    memorySample = (binornd(1,cueLevel)*2-1) + normrnd(0,1);
                end
                visualSample = (binornd(1,coherenceLevel)*2-1) + normrnd(0,1);
            else
                if ~mod(t,memoryThinning)
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
            
            if ~mod(t,memoryThinning)
                if memorySample > 0 
                    memSampsTargetA = 1 + memSampsTargetA;
                else
                    memSampsTargetB = 1 + memSampsTargetB;
                end
            end

            % update memory probability
            if order ==  2
                if ~mod(t,memoryThinning)
                    alphaMem = alphaMem + memSampsTargetA;
                    betaMem = betaMem + memSampsTargetB;
                    memoryAlphas(subj, trial, nSampMemory+t) = alphaMem;
                    memoryBetas(subj, trial, nSampMemory+t) = betaMem;
                    memoryPrecision = 1/betaVar(alphaMem, betaMem);
                end
            % and precision-weighted drift rate
                alphaVis = alphaVis + framesTargetA;
                betaVis = betaVis + framesTargetB;
                visualAlphas(subj, trial, t) = alphaVis;
                visualBetas(subj, trial, t) = betaVis;
                visualPrecision = 1/betaVar(alphaVis, betaVis);

            else % first-order analogs
                 % memory probability & precision
                 if ~mod(t,memoryThinning)
                    memProbTargetA = memSampsTargetA / (memSampsTargetA + memSampsTargetB);
                    memoryPrecision = 1/computeEntropy(memProbTargetA);
                 end
                 
                 % vision probability & precision
                 ghostProbTargetA = mean(ghostSamples(subj, trial, :) > 0);
                 probTargetA = (ghostProbTargetA + framesTargetA) / (v + framesTargetA + framesTargetB);
                 if t==1
                    visualPrecision = 1/computeEntropy(ghostProbTargetA);
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

             % additively accumulate precision-weighted evidence
             if t==1
                fullEvidence(subj, trial, t) = randn(1)+visualStartingPoint + memorySample*memoryDriftRate + visualSample*visualDriftRate;
             else
                 fullEvidence(subj, trial, t) = fullEvidence(subj, trial, t-1) + ...
                     memorySample*memoryDriftRate + visualSample*visualDriftRate;
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

%% save results
if saveFiles
    if order == 1
        %outfile = sprintf('results/first-order/workspace_files/%.2fcue_%.2fantCoh_%.2fcoh_%icong_%ithresh_%ithin_%iv_%imemSamp.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent, threshold, memoryThinning, v, nSampMemory);
         outfile = sprintf('ccn_submission/poster/supplement/FO_summary/%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    else
        %outfile = sprintf('results/second-order/workspace_files/%.2fcue_%.2fantCoh_%.2fcoh_%icong_%ithresh_%ithin_%iv_%imemSamp.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent, threshold, memoryThinning, v, nSampMemory);
         outfile = sprintf('ccn_submission/poster/supplement/SO_%.2fcue_%.2fantCoh_%.2fcoh_%icong.mat', cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    end
    save(outfile)

% write behavior csv for first order only
    if order == 1
        order = repmat(order, [nSub 1]);
        cueLevel = repmat(cueLevel, [nSub 1]);
        anticipatedCoherence = repmat(anticipatedCoherence, [nSub 1]);
        coherenceLevel = repmat(coherenceLevel, [nSub 1]);
        congruent = repmat(congruent, [nSub 1]);
        threshold = repmat(threshold, [nSub 1]);
        behav_csv = table(choices, RT, order, cueLevel, anticipatedCoherence, coherenceLevel, congruent, threshold);
        filename = ['ccn_submission/poster/supplement/FO_summary/' char(extractBetween(outfile, 'summary/', '.mat')) '.csv'];
%     if order == 1
%         %filename = ['results/first-order/behav_csv/' char(extractBetween(outfile, 'files/', '.mat')) '.csv'];
%     else
%         %filename = ['results/second-order/behav_csv/' char(extractBetween(outfile, 'files/', '.mat')) '.csv'];
%         filename = ['ccn_submission/poster/supplement/SO_' char(extractBetween(outfile, 'ccn', '.mat')) '.csv'];
%     end
        writetable(behav_csv, filename);
    end

    % write prior parameter csv
    if order == 2
        memoryAlphas = squeeze(memoryAlphas)';
        memoryBetas = squeeze(memoryBetas)';
        visualAlphas = [NaN(nSampMemory, nTrial); squeeze(visualAlphas)'];
        visualBetas = [NaN(nSampMemory, nTrial); squeeze(visualBetas)'];
        RT = repmat(RT, [110 1]);
        order = repmat(order, [nSampVisual + nSampMemory 1]);
        cueLevel = repmat(cueLevel, [nSampVisual + nSampMemory 1]);
        anticipatedCoherence = repmat(anticipatedCoherence, [nSampVisual + nSampMemory 1]);
        coherenceLevel = repmat(coherenceLevel, [nSampVisual + nSampMemory 1]);
        congruent = repmat(congruent, [nSampVisual + nSampMemory 1]);
        prior_params = table(memoryAlphas, memoryBetas, visualAlphas, visualBetas, RT, order, cueLevel, anticipatedCoherence, coherenceLevel, congruent);
        filename = ['ccn_submission/poster/supplement/priorParams_' char(extractBetween(outfile, 'SO_', '.mat')) '.csv'];
        writetable(prior_params, filename)
    end

end
%% Plot results
% trial traces
if plotResults
    if saveFiles
            cueLevel = unique(cueLevel);
            anticipatedCoherence = unique(anticipatedCoherence);
            coherenceLevel = unique(coherenceLevel);
            threshold = unique(threshold);
            order = unique(order);
    end
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
    if saveFiles
        figpath = ['ccn_submission/poster/supplement/' char(extractBetween(outfile, 'supplement/', '.mat'))];
        %outfig = char(regexp(outfile, 'b.*cong', 'match'));
        saveas(gcf, [figpath '_traces.png']);
    end
     
    % drifts
%     fig=figure; 
%     for b=1:10
%         subplot(2, 5, b)
%         hold on;
%         plot(squeeze(memoryDriftRates(1,b,:)), 'LineWidth',1.5)
%         plot([addX;squeeze(visualDriftRates(1, b, :))], 'LineWidth',1.5)
%         plot([1,140],[0,0],'k')
%         xline(nSampMemory, 'k--')
%         string2 = sprintf('trial %i', b);
%         title(string2);
%     end
%     h=legend({'memory','visual'},'Orientation','horizontal');
%     set(h, 'Position', [0.65 0.46 0.25 0.025]);
%     if order==1
%         sgtitle(sprintf(['first-order model\n' string]));
%     else
%         sgtitle(sprintf(['second-order model\n' string]));
%     end
%     plots=axes(fig, 'visible', 'off');
%     plots.XLabel.Visible='on';
%     plots.YLabel.Visible='on';
%     plots.Title.Visible='on';
%     xlabel(plots, 'time (a.u.)');
%     ylabel(plots, 'drift rate (a.u.)');
%     if saveFiles
%         saveas(gcf, [figpath '_drifts.png']);
%     end
end
end










