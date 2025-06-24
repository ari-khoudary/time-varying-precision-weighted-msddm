%% doSampling_summaryPlots
% requires the gramm toolbox for matlab bc i'm a ggplotter to my core

%% tidy data
if ~exist('summaryTable', 'var')
    tidyData;
end

%% raw & forced choice accuracy, proportion of no response trials, mean RT
clear f
% raw choice accuracy
f(1,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_rawChoice, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_rawChoice, 'ymax', summaryTable.upperCI_rawChoice);
f(1,1).set_names('x', 'coherence', 'y', 'propCorrect', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
f(1,1).set_title('raw choice (+/- sem)');
f(1,1).facet_grid(num2cell(num2str(summaryTable.memoryThinning), 2), num2cell(num2str(summaryTable.threshold), 2));
f(1,1).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(1,1).geom_point();
f(1,1).geom_interval('geom', 'errorbar');
f(1,1).geom_line();

% forced choice accuracy
f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_forcedChoice, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_forcedChoice, 'ymax', summaryTable.upperCI_forcedChoice);
f(1,2).set_names('x', 'coherence', 'y', 'propCorrect', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
f(1,2).set_title('forced choice (sign of accumulator) +/- sem');
f(1,2).facet_grid(num2cell(num2str(summaryTable.memoryThinning), 2), num2cell(num2str(summaryTable.threshold), 2));
f(1,2).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(1,2).geom_point();
f(1,2).geom_interval('geom', 'errorbar');
f(1,2).geom_line();
f(1,2).no_legend();


% proportion not hitting threshold
f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.nanRT,...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(2,1).set_names('x', 'coherence', 'y', 'count', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
f(2,1).set_title('number of trials with NaN choice');
f(2,1).facet_grid(num2cell(num2str(summaryTable.memoryThinning), 2), num2cell(num2str(summaryTable.threshold), 2));
f(2,1).geom_point();
f(2,1).geom_line();
f(2,1).no_legend();

% RT summary
f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_RT, ...
    'color', summaryTable.cue, 'linestyle', summaryTable.congruent,...
            'ymin', summaryTable.lowerCI_RT, 'ymax', summaryTable.upperCI_RT);
f(2,2).set_names('x', 'coherence', 'y', 'mean RT', ...
    'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
f(2,2).set_title('RT +/- sem; all trials');
f(2,2).facet_grid(num2cell(num2str(summaryTable.memoryThinning), 2), num2cell(num2str(summaryTable.threshold), 2));
f(2,2).geom_abline('slope', 0, 'intercept', 0.5, 'style', ':');
f(2,2).geom_abline('slope', 0, 'intercept', max(dataTable.nFrames), 'style', '-');
f(2,2).geom_point();
f(2,2).geom_interval('geom', 'errorbar');
f(2,2).geom_line();
f(2,2).no_legend();

f.set_line_options('style', {'-', '-.', '-'});
f.set_title(['n=' num2str(nSub) '; ' num2str(allData(1).nTrial) ' trials/cell/subj']);
figure('Name', 'Effects of threshold on accuracy & RT');
f.draw();

%% RT distributions for forced choice

% make vector of thinning levels
thinLevels = unique(dataTable.memoryThinning);
nFrames = dataTable{1,'nFrames'};
nBins = 15;

clear f
for i = 1:length(thinLevels)
    % draw distributions
    f(1,i) = gramm('x', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, 'subset', dataTable.memoryThinning==thinLevels(i) & dataTable.threshold<maxThresh);
    f(1,i).set_names('x', 'RT', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'c');
    f(1,i).set_title(['thinning= ' num2str(thinLevels(i))]);
    f(1,i).facet_grid(dataTable.coherence, dataTable.threshold, 'scale', 'independent');
    f(1,i).stat_bin('nbins', nBins, 'geom', 'line');

    % draw medians
    f(2,i) = gramm('x', dataTable.coherence, 'y', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
        'subset', dataTable.memoryThinning==thinLevels(i) & dataTable.threshold<maxThresh);
    f(2,i).set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh');
    f(2,i).set_title('median RT');
    f(2,i).facet_grid([], dataTable.threshold, 'scale', 'independent');
    f(2,i).stat_summary('type', '95percentile', 'geom', {'point', 'line'});
    f(2,i).no_legend();

     % draw means
    f(3,i) = gramm('x', dataTable.coherence, 'y', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
        'subset', dataTable.memoryThinning==thinLevels(i) & dataTable.threshold<maxThresh);
    f(3,i).set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh');
    f(3,i).set_title('mean RT');
    f(3,i).facet_grid([], dataTable.threshold, 'scale', 'independent');
    f(3,i).stat_summary('type', 'std', 'geom', {'point', 'line'});
    f(3,i).no_legend();

    if i<length(thinLevels)
        f(1,i).no_legend();
    end
end

f.set_line_options('style', {'-', '-.', '-'});
f.set_title('Forced choice RTs')
f.set_text_options('base_size', 8, 'legend_scaling', 1, 'legend_title_scaling', 1, 'facet_scaling', 1.2, 'title_scaling', 1.25);
figure('Name', 'forced choice RT distributions');
f.draw()

%% RT distributions for raw choice

clear f
for i = 1:length(thinLevels)
    % draw distributions
    f(1,i) = gramm('x', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, 'subset', dataTable.memoryThinning==thinLevels(i) & ~isnan(dataTable.rawChoice));
    f(1,i).set_names('x', 'RT', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'c');
    f(1,i).set_title(['thinning= ' num2str(thinLevels(i))]);
    f(1,i).facet_grid(dataTable.coherence, dataTable.threshold, 'scale', 'independent');
    f(1,i).stat_bin('nbins', nBins, 'geom', 'line');

    % draw medians
    f(2,i) = gramm('x', dataTable.coherence, 'y', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
        'subset', dataTable.memoryThinning==thinLevels(i) & ~isnan(dataTable.rawChoice));
    f(2,i).set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh');
    f(2,i).set_title('median RT');
    f(2,i).facet_grid([], dataTable.threshold, 'scale', 'independent');
    f(2,i).stat_summary('type', '95percentile', 'geom', {'point', 'line'});
    f(2,i).no_legend();

    % draw means
    f(3,i) = gramm('x', dataTable.coherence, 'y', dataTable.RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
        'subset', dataTable.memoryThinning==thinLevels(i) & ~isnan(dataTable.rawChoice));
    f(3,i).set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh');
    f(3,i).set_title('mean RT');
    f(3,i).facet_grid([], dataTable.threshold, 'scale', 'independent');
    f(3,i).stat_summary('type', 'std', 'geom', {'point', 'line'});
    f(3,i).no_legend();

    % only put a legend for the last distribution plot
    if i<length(thinLevels)
        f(1,i).no_legend();
    end
end

f.set_line_options('style', {'-', '-.', '-'});
f.set_title('Raw choice RTs')
f.set_text_options('base_size', 8, 'legend_scaling', 1, 'legend_title_scaling', 1, 'facet_scaling', 1.2, 'title_scaling', 1.25);
figure('Name', 'raw choice RT distributions');
f.draw()

%% starting point distributions & medians

% plot distributions 
clear f
f(1,1) = gramm('x', dataTable.dummyStart1, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
f(1,1).facet_grid(dataTable.coherence, dataTable.memoryThinning, 'scale', 'independent');
f(1,1).set_names('x', 'dv value', 'y', 'density', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'coh');
f(1,1).set_title('starting point (note axis scale change across plots)');
f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
f(1,1).stat_bin('nbins', 30, 'geom', 'line');

% plot medians
f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.dummyMedian_startPoint1, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(1,2).facet_grid([], summaryTable.memoryThinning);
f(1,2).set_names('x', 'coherence', 'y', 'median start point', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thinning');
f(1,2).set_title('median averaged across threshold (+/- sem)');
f(1,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
f(1,2).stat_summary('type', 'sem', 'geom', {'point', 'errorbar', 'line'});

f.set_line_options('style', {'-', '-.', '-'});
figure('Name', 'DV starting point summaries');
f.draw()

%% drift rate distributions & means

figNames = {'slope of DV during first noise period', 'slope of DV during first signal period', 'slope of DV during second noise period', 'slope of DV during second signal period'};

for i = 1:length(figNames)
    clear f
    % initialize a different figure for each epoch of the experiment
    if i==1 % first noise period
        % distribution across thresholds
        f(1,1) = gramm('x', dataTable.dummySlope1, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.dummyMedian_dvSlope1, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
    elseif i==2 % first signal period
        f(1,1) = gramm('x', dataTable.dummySlope2, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.dummyMedian_dvSlope2, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
    elseif i==3 % second noise period
        f(1,1) = gramm('x', dataTable.dummySlope3, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.dummyMedian_dvSlope3, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
    else % second signal period
        f(1,1) = gramm('x', dataTable.dummySlope4, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.dummyMedian_dvSlope4, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
    end


    % draw distributions
    f(1,1).set_title('slope values across thresholds (note scale shifts)');
    f(1,1).set_names('x', 'dv slope', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'coh');
    f(1,1).facet_grid(dataTable.coherence, dataTable.memoryThinning, 'scale', 'independent');
    f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
    f(1,1).stat_bin('nbins', 15, 'geom', 'line');
    f(1,1).no_legend();

    % draw medians
    f(1,2).set_title('median slope');
    f(1,2).set_names('x', 'coh', 'y', 'median slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'threshold');
    f(1,2).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
    f(1,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
    f(1,2).geom_point();
    f(1,2).geom_line();

    f.set_title(figNames{i});
    f.set_line_options('style', {'-', '-.', '-'});
    f.set_text_options('legend_scaling', 0.8, 'legend_title_scaling', 0.8);
    figure('Name', 'Effect of conditions on slope of DV');
    f.draw();
end

%% plot evidence-relative RT distributions for all trials
clear f
% first noise
f(1,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_noise1RT, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(1,1).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
f(1,1).set_title('first noise period');
f(1,1).geom_point();
f(1,1).geom_line();
% first signal
f(1,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_signal1RT, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(1,2).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
f(1,2).set_title('first signal period');
f(1,2).geom_point();
f(1,2).geom_line();
% second noise
f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_noise2RT, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(2,1).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
f(2,1).set_title('second noise period');
f(2,1).geom_point();
f(2,1).geom_line();
% second signal
f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_signal2RT, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent);
f(2,2).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
f(2,2).set_title('second signal period');
f(2,2).geom_point();
f(2,2).geom_line();

f.set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'row', 'thresh', 'column', 'thin');
f.set_line_options('style', {'-', '-.', '-'});
f.set_text_options('base_size', 8, 'legend_scaling', 1, 'legend_title_scaling', 1, 'facet_scaling', 1.2, 'title_scaling', 1.25);
f.set_title('proportion of RTs in each trial epoch (forced choice); note y-axis difference across plots')
figure('Name', 'evidence-relative forced choice responses');
f.draw();

%% plot evidence-relative RT distributions for raw choice trials
clear f
% first noise
f(1,1) = gramm('x', dataTable.coherence, 'y', dataTable.noise1RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
    'subset', ~isnan(dataTable.rawChoice));
f(1,1).facet_grid(dataTable.threshold, dataTable.memoryThinning);
f(1,1).set_title('first noise period');
f(1,1).stat_summary('type', 'std', 'geom', {'point', 'line'});
% first signal
f(1,2) = gramm('x', dataTable.coherence, 'y', dataTable.signal1RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
    'subset', ~isnan(dataTable.rawChoice));
f(1,2).facet_grid(dataTable.threshold, dataTable.memoryThinning);
f(1,2).set_title('first signal period');
f(1,2).stat_summary('type', 'std', 'geom', {'point', 'line'});
% second noise
f(2,1) = gramm('x', dataTable.coherence, 'y', dataTable.noise2RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
    'subset', ~isnan(dataTable.rawChoice));
f(2,1).facet_grid(dataTable.threshold, dataTable.memoryThinning);
f(2,1).set_title('second noise period');
f(2,1).stat_summary('type', 'std', 'geom', {'point', 'line'});
% second signal
f(2,2) = gramm('x', dataTable.coherence, 'y', dataTable.signal2RT, 'color', dataTable.cue, 'linestyle', dataTable.congruent, ...
    'subset', ~isnan(dataTable.rawChoice));
f(2,2).facet_grid(dataTable.threshold, dataTable.memoryThinning);
f(2,2).set_title('second signal period');
f(2,2).stat_summary('type', 'std', 'geom', {'point', 'line'});

f.set_names('x', 'coh', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'row', 'thresh', 'column', 'thin');
f.set_line_options('style', {'-', '-.', '-'});
f.set_text_options('base_size', 8, 'legend_scaling', 1, 'legend_title_scaling', 1, 'facet_scaling', 1.2, 'title_scaling', 1.25);
f.set_title('proportion of RTs in each trial epoch (raw choice); note y-axis difference across plots')
figure('Name', 'evidence-relative raw choice responses');
f.draw();

%% plot relative RT as stacked histograms
clear f
f(1,1) = gramm('x', dataTable.coherence, 'color', dataTable.relativeRT, 'subset', ~isnan(dataTable.rawChoice));
f(1,1).facet_grid(dataTable.cue, dataTable.congruent, 'scale', 'independent')
f(1,1).stat_bin('geom', 'stacked_bar', 'width', 1)
f.set_names('x', 'coherence', 'y', 'count', 'column', '', 'row', 'cue');
f.coord_flip();
figure;
f.draw()

%% plot effects of time on starting point offset
clear f
f(1,1) = gramm('x', dataTable.dummyStart1, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
f(1,1).facet_grid(dataTable.memoryThinning, dataTable.noise1Quartile, 'scale', 'free');
f(1,1).set_names('x', 'dv value', 'y', 'density', 'color', 'cue', 'linestyle', 'congruent', 'column', 'noise1 quartile', 'row', 'thin');
f(1,1).set_title('starting point (note axis scale change across plots)');
f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
f(1,1).stat_bin('nbins', 30, 'geom', 'line');

% plot medians
f(1,2) = gramm('x', tv_summaryTable.noise1Quartile, 'y', tv_summaryTable.dummyMedian_startPoint1, 'color', tv_summaryTable.cue, 'linestyle', tv_summaryTable.congruent);
f(1,2).facet_grid([], tv_summaryTable.memoryThinning);
f(1,2).set_names('x', 'noise 1 quartile', 'y', 'median start point', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thinning');
f(1,2).set_title('median averaged across threshold (+/- sem)');
f(1,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
f(1,2).stat_summary('type', 'sem', 'geom', {'point', 'errorbar', 'line'});

f.set_line_options('style', {'-', '-.', '-'});
figure('Name', 'starting point offset as a function of time');
f.draw()


%% plot effects of time on dv slope

tv_figNames = {'slope of DV during first noise period', 'slope of DV during second noise period'};
for i = 1:length(tv_figNames)
    clear f
    % initialize a different figure for each epoch of the experiment
    if i==1 % first noise period
        f(1,1) = gramm('x', dataTable.dummySlope1, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', tv_summaryTable.noise1Quartile, 'y', tv_summaryTable.dummyMean_dvSlope1, 'color', tv_summaryTable.cue, 'linestyle', tv_summaryTable.congruent);
    else % second noise period
        f(1,1) = gramm('x', dataTable.dummySlope3, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
        f(1,2) = gramm('x', tv_summaryTable.noise2Quartile, 'y', tv_summaryTable.dummyMean_dvSlope3, 'color', tv_summaryTable.cue, 'linestyle', tv_summaryTable.congruent);
    end

    % draw distributions
    f(1,1).set_title('slope values across thresholds (note scale shifts)');
    if i==1
        f(1,1).set_names('x', 'dv slope', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'noise1 quartile', 'row', 'thin');
        f(1,1).facet_grid(dataTable.memoryThinning, dataTable.noise1Quartile, 'scale', 'free');
    else
        f(1,1).set_names('x', 'dv slope', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'noise2 quartile', 'row', 'thin');
        f(1,1).facet_grid(dataTable.memoryThinning, dataTable.noise2Quartile, 'scale', 'free');
    end
    f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
    f(1,1).stat_bin('nbins', 15, 'geom', 'line');
    f(1,1).no_legend();

    % draw medians
    f(1,2).set_title('median slope');
    if i==1
        f(1,2).set_names('x', 'noise 1 duration', 'y', 'mean slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
        f(1,2).facet_grid(tv_summaryTable.memoryThinning, tv_summaryTable.threshold, 'scale', 'free_y');
    else
        f(1,2).set_names('x', 'coh', 'y', 'median slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thresh', 'row', 'thin');
        f(1,2).facet_grid(tv_summaryTable.memoryThinning, tv_summaryTable.threshold, 'scale', 'free_y');
    end
    f(1,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
    % average across coherence levels
    f(1,2).stat_summary('type', 'ci', 'geom', {'point', 'line'});

    f.set_title(tv_figNames{i});
    f.set_line_options('style', {'-', '-.', '-'});
    f.set_text_options('legend_scaling', 0.8, 'legend_title_scaling', 0.8);
    figure('Name', 'Effect of time on slope of DV');
    f.draw();
end


%% effect of signal 1 duration on slope

signal1quartiles = quantile(dataTable.signal1Frames, 4);
signal2quartiles = quantile(dataTable.signal2Frames, 4);
dataTable.signal1Quartile = zeros(height(dataTable), 1);
dataTable.signal2Quartile = dataTable.signal1Quartile;

for i=1:4
    if i==1
        dataTable.signal1Quartile(dataTable.signal1Frames <= signal1quartiles(i)) = i;
        dataTable.signal2Quartile(dataTable.signal2Frames <= signal2quartiles(i)) = i;
    elseif i<4
        dataTable.signal1Quartile(dataTable.signal1Frames <= signal1quartiles(i) & dataTable.signal1Frames > signal1quartiles(i-1)) = i;
        dataTable.signal2Quartile(dataTable.signal2Frames <= signal2quartiles(i) & dataTable.signal2Frames > signal2quartiles(i-1)) = i;
    else
        dataTable.signal1Quartile(dataTable.signal1Frames >= signal1quartiles(i)) = i;
        dataTable.signal2Quartile(dataTable.signal2Frames >= signal2quartiles(i)) = i;
    end
end

dataTable.signal1Quartile(dataTable.signal1Frames <= signal1quartiles(1)) = 1;
dataTable.signal1Quartile(dataTable.signal1Frames > signal1quartiles(1) & dataTable.signal1Frames < signal1quartiles(3)) = 2;
dataTable.signal1Quartile(dataTable.signal1Frames < signal1quartiles(4) & dataTable.signal1Frames >= signal1quartiles(3)) = 3;
dataTable.signal1Quartile(dataTable.signal1Frames >= signal1quartiles(4)) = 4;

clear f
f(1,1) = gramm('x', dataTable.dummySlope2, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
f(1,1).set_title(['slope values for thresh=' num2str(maxThresh)]);
f(1,1).set_names('x', 'dv slope', 'y', '', 'color', 'cue', 'linestyle', 'congruent', 'column', 'signal 1 quartile', 'row', 'thin');
f(1,1).facet_grid(dataTable.memoryThinning, dataTable.signal1Quartile, 'scale', 'free');
f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
f(1,1).stat_bin('nbins', 15, 'geom', 'line');
f(1,1).no_legend();

f(1,2) = gramm('x', tv_summaryTable.signal1Quartile, 'y', tv_summaryTable.dummyMean_dvSlope2, 'color', tv_summaryTable.cue, 'linestyle', tv_summaryTable.congruent);
f(1,2).set_title('median slope');
f(1,2).set_names('x', 'signal 1 quartile', 'y', 'mean slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'coh', 'row', 'thin');
f(1,2).facet_grid(tv_summaryTable.memoryThinning, tv_summaryTable.coherence, 'scale', 'free_y');
f(1,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
% average across coherence levels
f(1,2).stat_summary('type', 'ci', 'geom', {'point', 'line'});

f.set_title('slope of dv during first signal period')
f.set_line_options('style', {'-', '-.', '-'});
f.set_text_options('legend_scaling', 0.8, 'legend_title_scaling', 0.8);
figure('Name', 'Effect of time on slope of DV');
f.draw();

%% plot effects of threshold & thining on dv slope (medians & means)
% 
% for i = 3
%     clear f
%     % initialize a different figure for each epoch of the experiment
%     if i==1 % first noise period
%         % distribution across thresholds
%         f(1,1) = gramm('x', dataTable.dvSlope1, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
%         % median +/- sem
%         f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.median_dvSlope1, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.median_dvSlope1 - summaryTable.fun1_dvSlope1, 'ymax', summaryTable.median_dvSlope1 + summaryTable.fun1_dvSlope1);
%         % mean +/- sem
%         f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_dvSlope1, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.mean_dvSlope1 - summaryTable.fun1_dvSlope1, 'ymax', summaryTable.mean_dvSlope1 + summaryTable.fun1_dvSlope1);
%         f(1,1).set_title('DV slope - first noise period (all thresholds)');
%     elseif i==2 % first signal period
%         % distribution across thresholds
%         f(1,1) = gramm('x', dataTable.dvSlope2, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
%         % median +/- sem
%         f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.median_dvSlope2, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.median_dvSlope2 - summaryTable.fun1_dvSlope2, 'ymax', summaryTable.median_dvSlope2 + summaryTable.fun1_dvSlope2);
%         % mean +/- sem
%         f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_dvSlope2, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.mean_dvSlope2 - summaryTable.fun1_dvSlope2, 'ymax', summaryTable.mean_dvSlope2 + summaryTable.fun1_dvSlope2);
%         f(1,1).set_title('DV slope - first signal period (all thresholds)');
%     elseif i==3 % second noise period
%         % distribution across thresholds
%         f(1,1) = gramm('x', dataTable.dvSlope3, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
%         % median +/- sem
%         f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.median_dvSlope3, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.median_dvSlope3 - summaryTable.fun1_dvSlope3, 'ymax', summaryTable.median_dvSlope3 + summaryTable.fun1_dvSlope3);
%         % mean +/- sem
%         f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_dvSlope3, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.mean_dvSlope3 - summaryTable.fun1_dvSlope3, 'ymax', summaryTable.mean_dvSlope3 + summaryTable.fun1_dvSlope3);
%         f(1,1).set_title('DV slope - second noise period (all thresholds)');
%     else % second signal period
%         % distribution
%         f(1,1) = gramm('x', dataTable.dvSlope4, 'color', dataTable.cue, 'linestyle', dataTable.congruent);
%         % median +/- sem
%         f(2,1) = gramm('x', summaryTable.coherence, 'y', summaryTable.median_dvSlope4, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.median_dvSlope4 - summaryTable.fun1_dvSlope4, 'ymax', summaryTable.median_dvSlope4 + summaryTable.fun1_dvSlope4);
%         % mean +/- sem
%         f(2,2) = gramm('x', summaryTable.coherence, 'y', summaryTable.mean_dvSlope4, 'color', summaryTable.cue, 'linestyle', summaryTable.congruent, ...
%             'ymin', summaryTable.mean_dvSlope4 - summaryTable.fun1_dvSlope4, 'ymax', summaryTable.mean_dvSlope4 + summaryTable.fun1_dvSlope4);
%         f(1,1).set_title('DV slope - second signal period (all thresholds)');
%     end
% 
%     % draw distributions
%     f(1,1).set_names('x', 'dv slope', 'y', 'density', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'coh');
%     f(1,1).facet_grid(dataTable.coherence, dataTable.memoryThinning);
%     f(1,1).geom_vline('xintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
%     f(1,1).stat_bin('nbins', 30, 'geom', 'line');
% 
%     % draw medians
%     f(2,1).set_names('x', 'coh', 'y', 'median slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'threshold');
%     f(2,1).set_title('median slope (+/- sem)');
%     f(2,1).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
%     f(2,1).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
%     f(2,1).geom_point();
%     f(2,1).geom_interval('geom', 'errorbar');
%     f(2,1).geom_line();
% 
%     % draw means
%     f(2,2).set_names('x', 'coh', 'y', 'mean slope', 'color', 'cue', 'linestyle', 'congruent', 'column', 'thin', 'row', 'threshold');
%     f(2,2).set_title('mean slope +/- sem');
%     f(2,2).facet_grid(summaryTable.threshold, summaryTable.memoryThinning);
%     f(2,2).geom_hline('yintercept', 0, 'style', 'k-', 'LineWidth', 0.5);
%     f(2,2).geom_point();
%     f(2,2).geom_interval('geom', 'errorbar');
%     f(2,2).geom_line();
% 
%     f.set_line_options('style', {'-', '-.', '-'})
%     figure('Name', 'Effect of conditions on slope of DV');
%     f.draw();
% end
% 
% 
