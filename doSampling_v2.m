function [choices, RT] = doSampling_v2(cueLevel, coherence, congruent, threshold, memoryThinning, ...
    saveFiles, plotResults)

% multi-stage, multi-source sequential sampling model. sample weights vary
% dynamically as the relative precision of each evidence stream, updated
% after each sample.

% INPUTS
% cueLevel [0:1]: sets the precision of memory evidence

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

% define simulation parameters
nSub = 1;
nTrial = 1000;
%coherence = 0.8;
visualThinning = 1;
trialDuration = 3000/2; % maximum duration of evidence stream divided by 2 to account for inter-frame noise

% create variables to store values
choices = zeros(nTrial, 1);
RT = zeros(nTrial, 1);
noise1_duration= zeros(nTrial, 1);
noise2_duration = zeros(nTrial, 1);
noise2_onset = zeros(nTrial, 1);
samples = zeros(trialDuration, 2, nTrial);
evidence = zeros(trialDuration, 3, nTrial);
driftRates = zeros(trialDuration, 2, nTrial);
precisions = driftRates;
counters = zeros(trialDuration, 4, nTrial);

% precompute entropy distribution
p = 0.05:.1:0.95;
entropy = zeros(length(p), 1);
for i = 1:length(p)
    entropy(i) = computeEntropy(p(i));
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
    for trial=1:nTrial
        % initialize trial structure
        noise1Duration = round(unifrnd(250/2, 1000/2));
        noise2Duration = round(unifrnd(250/2, 1000/2));
        noise2Onset = round(unifrnd(noise1Duration + trialDuration*0.1, trialDuration-noise2Duration));
        signal1Duration = noise2Onset - noise1Duration;
        signal2Duration = trialDuration - (noise2Duration+noise2Onset);
        signal2Onset = noise1Duration + signal1Duration + noise2Duration + 1;

        % save noise durations
        noise1_duration(trial) = noise1Duration;
        noise2_onset(trial) = noise2Onset;
        noise2_duration(trial) = noise2Duration;

        % initialize counters
        alphaMem = 1;
        betaMem = 1;
        alphaVis = 1;
        betaVis = 1;

        % generate & accumulate memory evidence samples
        samples(:, memory_idx, trial) = (binornd(1,cueLevel, [trialDuration,1])*2-1) + normrnd(0,1, [trialDuration,1]);

        % generate visual evidence samples during signal periods only
        if congruent
            samples(noise1Duration+1:noise2Onset, visual_idx, trial) = (binornd(1,coherence, [signal1Duration,1])*2-1) + normrnd(0,1,[signal1Duration,1]);
            samples(signal2Onset:trialDuration, visual_idx, trial) = (binornd(1,coherence, [signal2Duration,1])*2-1) + normrnd(0,1,[signal2Duration,1]);
        else
            samples(noise1Duration+1:noise2Onset, visual_idx, trial) = -(binornd(1,coherence, [signal1Duration,1])*2-1) + normrnd(0,1,[signal1Duration,1]);
            samples(signal2Onset:trialDuration, visual_idx, trial) = -(binornd(1,coherence, [signal2Duration,1])*2-1) + normrnd(0,1,[signal2Duration,1]);
        end
        % populate with noise otherwise
        samples(1:noise1Duration, visual_idx, trial) = normrnd(0,1,[noise1Duration,1]);
        samples(noise2Onset:noise2Onset+noise2Duration, visual_idx, trial) = normrnd(0,1,[noise2Duration+1,1]);

        % compute time-varying drift rate
        for t=1:trialDuration
            if ~mod(t, visualThinning)
                % update visual counters
                if samples(t, visual_idx, trial) > 0
                    alphaVis = alphaVis + 1;
                    counters(t, alphaVis_idx, trial) = alphaVis;
                else
                    betaVis = betaVis + 1;
                    counters(t, betaVis_idx, trial) = betaVis;
                end

                % compute normalized beta pdf
                visualPDF = betapdf(p, alphaVis, betaVis) ./ sum(betapdf(p, alphaVis, betaVis));
                % compute precision as inverse belief-weighted entropy
                visualPrecision = 1/sum(visualPDF * entropy);
                precisions(t, visual_idx, trial) = visualPrecision;
            end

            if ~mod(t, memoryThinning)
                % update memory counters
                if samples(t, memory_idx, trial) > 0
                    alphaMem = alphaMem + 1;
                    counters(t, alphaMem_idx, trial) = alphaMem;
                else
                    betaMem = betaMem + 1;
                    counters(t, betaMem_idx, trial) = betaMem;
                end                

                % compute normalized beta pdf
                memoryPDF = betapdf(p, alphaMem, betaMem) ./ sum(betapdf(p, alphaMem, betaMem));
                % compute precision as inverse belief-weighted entropy
                memoryPrecision = 1/sum(memoryPDF * entropy);
                precisions(t, memory_idx, trial) = memoryPrecision;
            
            elseif t < memoryThinning
                    memoryPrecision = 0;
            else
                    prevPrecision = find(precisions(:, memory_idx, trial), 1, 'last');
                    memoryPrecision = precisions(prevPrecision, memory_idx, trial);
            end

            % compute relative precision evidence weights
            visualDriftRate = visualPrecision / (visualPrecision + memoryPrecision);
            memoryDriftRate = memoryPrecision / (visualPrecision + memoryPrecision);
            driftRates(t, visual_idx, trial) = visualDriftRate;
            driftRates(t, memory_idx, trial) = memoryDriftRate;

            % compute time-varying relative precision-weighted decision variable
            memorySample = samples(t, memory_idx, trial);
            visualSample = samples(t, visual_idx, trial);
            
            if t == 1
                evidence(t, visual_idx, trial) = visualSample * visualDriftRate;
                if ~mod(t, memoryThinning)
                    evidence(t, memory_idx, trial) = memorySample * memoryDriftRate;
                    evidence(t, dv_idx, trial) = memorySample*memoryDriftRate + visualSample * visualDriftRate;
                else
                    evidence(t, dv_idx, trial) = evidence(t, dv_idx, trial) + visualSample*visualDriftRate;
                end

            else
                evidence(t, visual_idx, trial) = evidence(t-1, visual_idx, trial) + visualSample * visualDriftRate;
                if ~mod(t, memoryThinning)
                    evidence(t, memory_idx, trial) = evidence(t-1, memory_idx, trial) + memorySample * memoryDriftRate;
                    evidence(t, dv_idx, trial) = evidence(t-1, dv_idx, trial) + memorySample*memoryDriftRate + visualSample * visualDriftRate;
                else
                    evidence(t, memory_idx, trial) = evidence(t-1, memory_idx, trial);
                    evidence(t, dv_idx, trial) = evidence(t-1, dv_idx, trial) + visualSample*visualDriftRate;
                end
            end
        end

        % find point at which evidence crosses threshold
        dv = evidence(:, dv_idx, trial);
        boundaryIdx = find(abs(dv)>threshold, 1);
        if isempty(boundaryIdx)
            boundaryIdx = trialDuration;
        end

        % store crossing point as RT
        RT(trial) = boundaryIdx;

        % populate choice matrix accordingly
        if congruent 
            if dv(boundaryIdx) > threshold
                choices(trial, 1) = 1;
            else
                choices(trial, 1) = 0;
            end
        else
            if dv(boundaryIdx) < -threshold
                choices(trial, 1) = 1;
            else
                choices(trial, 1) = 0;
            end
        end
    end
end

%% save results
if saveFiles
    % save workspace
    outfile = sprintf('results_v2/workspace_files/equalN1_%.2fcue_%icong_%.2fcoh_%ithresh.mat', cueLevel, congruent, coherence, threshold);
    save(outfile)

    % write CSVs
    behav_csv = table(choices, RT, noise1_duration, noise2_onset, noise2_duration);
    csv_name = sprintf('results_v2/behav_csv/equalN1_%.2fcue_%icong_%.2fcoh_%ithresh.csv', cueLevel, congruent, coherence, threshold);
    writetable(behav_csv, csv_name);

    %     % write prior parameter csv
    %     if order == 2
    %         memoryAlphas = squeeze(memoryAlphas)';
    %         memoryBetas = squeeze(memoryBetas)';
    %         visualAlphas = [NaN(nSampMemory, nTrial); squeeze(visualAlphas)'];
    %         visualBetas = [NaN(nSampMemory, nTrial); squeeze(visualBetas)'];
    %         RT = repmat(RT, [110 1]);
    %         order = repmat(order, [trialDuration + nSampMemory 1]);
    %         cueLevel = repmat(cueLevel, [trialDuration + nSampMemory 1]);
    %         anticipatedCoherence = repmat(anticipatedCoherence, [trialDuration + nSampMemory 1]);
    %         coherenceLevel = repmat(coherenceLevel, [trialDuration + nSampMemory 1]);
    %         congruent = repmat(congruent, [trialDuration + nSampMemory 1]);
    %         prior_params = table(memoryAlphas, memoryBetas, visualAlphas, visualBetas, RT, order, cueLevel, anticipatedCoherence, coherenceLevel, congruent);
    %         filename = ['ccn_submission/poster/supplement/priorParams_' char(extractBetween(outfile, 'SO_', '.mat')) '.csv'];
    %         writetable(prior_params, filename)
    %     end
end

%% Plot results
% trial traces
if plotResults
    fig=figure;
    for a=1:10
        subplot(2, 5, a)
        hold on;
        plot(evidence(:, memory_idx, a), 'LineWidth',1.5)
        plot(evidence(:, visual_idx, a), 'LineWidth',1.5)
        plot(evidence(:, 3, a), 'LineWidth',1.5)
        xline(noise1_duration(a))
        xline([noise2_onset(a), noise2_onset(a)+noise2_duration(a)], 'k--')
        yline([threshold -threshold])
        %ylim([-threshold-500, threshold+500])
        string2 = sprintf('trial %i', a);
        title(string2);
    end

    % add title
    if congruent
        string = sprintf('%.2f cue, %.2f coherece, %i thinning, congruent', cueLevel, coherence, memoryThinning);
    else
        string = sprintf('%.2f cue, %.2f coherence, %i thinning, incongruent', cueLevel, coherence, memoryThinning);
    end
    sgtitle(sprintf(['combo model\n' string]));

    % make pretty
    h=legend({'memory','visual','combined'},'FontSize',8, 'Orientation', 'horizontal');
    set(h, 'Position', [0.65 0.46 0.25 0.025]);
    plots=axes(fig, 'visible', 'off');
    plots.XLabel.Visible='on';
    plots.YLabel.Visible='on';
    plots.Title.Visible='on';
    xlabel(plots, 'time (a.u.)');
    ylabel(plots, 'evidence (a.u.)');
    if saveFiles
        figpath = sprintf('results_v2/trace_figs/%.2fcue_%icong_%.2fcoh_%ithin.png', cueLevel, congruent, coherence, memoryThinning);
        saveas(gcf, figpath);
    end

    % drifts
%     fig=figure;
%     for b=1:10
%         subplot(2, 5, b)
%         hold on;
%         plot(driftRates(:, memory_idx, b), 'LineWidth',1.5)
%         plot(driftRates(:, visual_idx, b), 'LineWidth',1.5)
%         xline(noise1_duration(b))
%         xline([noise2_onset(b), noise2_onset(b)+noise2_duration(b)], 'k--')
%         string2 = sprintf('trial %i', b);
%         title(string2);
%     end
% 
%     % add title
%     if congruent
%         string = sprintf('%.2f cue, %.2f coherece, %i thinning, congruent', cueLevel, coherence, memoryThinning);
%     else
%         string = sprintf('%.2f cue, %.2f coherece, %i thinning, incongruent', cueLevel, coherence, memoryThinning);
%     end
%     sgtitle(sprintf(['combo model\n' string]));
% 
%     % make pretty
%     h=legend({'memory','visual'},'Orientation','horizontal');
%     set(h, 'Position', [0.65 0.46 0.25 0.025]);
%     plots=axes(fig, 'visible', 'off');
%     plots.XLabel.Visible='on';
%     plots.YLabel.Visible='on';
%     plots.Title.Visible='on';
%     xlabel(plots, 'time (a.u.)');
%     ylabel(plots, 'drift rate (a.u.)');
%     if saveFiles
%         figpath = sprintf('results_v2/drift_figs/%.2fcue_%icong_%.2fcoh_%ithin.png', cueLevel, congruent, coherence, memoryThinning);
%         saveas(gcf, figpath);
%     end
end
end










