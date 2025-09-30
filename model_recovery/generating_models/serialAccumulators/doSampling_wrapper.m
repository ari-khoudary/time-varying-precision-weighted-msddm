%% specify simulation settings
clear
debug = 0;
cue = [0.5 0.65 0.8];
threshold_coherence = readmatrix("../calibratedCoherences.csv");
threshold_coherence = threshold_coherence(:, 1:2);
memoryThinning = 4:4:16;
visionThinning = 1;
vizPresentationRate = 1/60;
if ~debug
    nSub = 50;
    nTrial = 1000; 
else
    nSub = 1;
    nTrial = 100;
end

% noise periods -- correspond to expt design as of 08-2025
noisePeriods = 1; % logical: do you want 2 noise periods?
noNoiseTrialDuration = 3; % in seconds: how long do you want trials to be if there are no noise periods?
% parameters of the noise periods
expLambda = 0.12; % parameter of the exponential defining the hazard rates 
minNoiseDuration = 0.75; % seconds
minSignalDuration = 0.5; % seconds

% half neutral trials (boolean)
halfNeutralTrials = 0;

% visual evidence noise
flickerAdditiveNoise = 1;  % logical; do you want to add noise to each sample of visual evidence?
flickerAdditiveNoiseValue = 'gaussian';  % string; what kind of noise do you want to add to each visual evidence sample? gaussian=zero-centered gaussian
flickerPadding = 1;  % logical; do you want to pad each signal frame with a noise frame?
flickerPaddingValue = 'zero'; % string; how do you want to model noise frames? options are zeros, zero-centered gaussian, more to come

% where do you want to save the results? (subdirectory of current dir)
timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd'));
outDir = sprintf('results/%s/', timestamp);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Create log directory
logDir = sprintf('%s/logs/', outDir);
if ~exist(logDir, 'dir')
    mkdir(logDir);
end

%% Get array job parameters
% Check if running as array job
if isempty(getenv('SLURM_ARRAY_TASK_ID'))
    % Running locally or single job - process all subjects
    startSub = 1;
    endSub = nSub;
    jobID = 1;
else
    % Running as array job
    arrayID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    jobID = arrayID;
    
    % Determine subjects per job
    subsPerJob = 10; % Process 10 subjects per array task
    startSub = (arrayID - 1) * subsPerJob + 1;
    endSub = min(arrayID * subsPerJob, nSub);
    
    fprintf('Array job %d processing subjects %d to %d\n', arrayID, startSub, endSub);
end

%% create config files
nCombo = length(threshold_coherence) * length(memoryThinning) * length(cue);
allConfigs = repmat({struct('myfield', {})}, 1, nCombo);

counter=0;
for a = 1:length(threshold_coherence)
    for b = 1:length(cue)
        for s = 1:length(memoryThinning)

            counter=counter+1;

            config.nTrial = nTrial;
            config.threshold = threshold_coherence(a, 1);
            config.coherence = threshold_coherence(a, 2);
            config.cue = cue(b);
            config.memoryThinning = memoryThinning(s);
            config.visionThinning = visionThinning;
            config.vizPresentationRate = vizPresentationRate;
            config.noisePeriods = noisePeriods;
            config.noNoiseTrialDuration = noNoiseTrialDuration;

            config.minNoiseDuration = minNoiseDuration;
            config.minSignalDuration = minSignalDuration;
            config.expLambda = expLambda;

            config.halfNeutralTrials = halfNeutralTrials;
            config.flickerAdditiveNoise = flickerAdditiveNoise;
            config.flickerAdditiveNoiseValue = flickerAdditiveNoiseValue;
            config.flickerNoisePadding = flickerPadding;
            config.flickerPaddingValue = flickerPaddingValue;
            config.outDir = outDir;
            config.logDir = logDir;

            allConfigs{counter} = config;
        end
    end
end

%% Initialize progress tracking
progressFile = sprintf('%s/progress_job_%d.mat', logDir, jobID);
completedFile = sprintf('%s/completed_subjects.txt', logDir);

% Check for existing progress
if exist(progressFile, 'file')
    load(progressFile, 'completedSubjects', 'completedConfigs');
    fprintf('Resuming from previous run. Found %d completed subjects.\n', length(completedSubjects));
else
    completedSubjects = [];
    completedConfigs = [];
end

%% Process subjects in this job
totalSubsThisJob = endSub - startSub + 1;
processedSubsThisJob = 0;

% start parallel pool
if ~debug && totalSubsThisJob > 1
    if isempty(gcp('nocreate'))
        parpool('local');
    end
end

%% Main processing loop
for subj = startSub:endSub
    % Check if this subject is already completed
    if ismember(subj, completedSubjects)
        fprintf('Subject %d already completed, skipping.\n', subj);
        continue;
    end
    
    fprintf('Processing subject %d of %d (job %d)\n', subj, nSub, jobID);
    
    % Process all configs for this subject
    subjectData = cell(length(allConfigs), 1);
    
    % Run all configurations for this subject in parallel
    parfor configIdx = 1:length(allConfigs)
        config = allConfigs{configIdx};
        config.subID = subj;
        subjectData{configIdx} = doSampling_serial(config);
    end
    
    % Combine all data for this subject
    combinedData = [];
    for configIdx = 1:length(allConfigs)
        if isempty(combinedData)
            combinedData = subjectData{configIdx};
        else
            combinedData = [combinedData; subjectData{configIdx}];
        end
    end
    
    % Write subject data immediately
    outfile = sprintf('%s/sub%d.csv', outDir, subj);
    writetable(combinedData, outfile);
    
    % Update progress tracking
    completedSubjects = [completedSubjects, subj];
    save(progressFile, 'completedSubjects', 'completedConfigs');
    
    % Log completion
    logEntry = sprintf('%s: Job %d completed subject %d\n', datetime('now'), jobID, subj);
    appendToLog(completedFile, logEntry);
    
    processedSubsThisJob = processedSubsThisJob + 1;
    fprintf('Completed subject %d (%d/%d in this job)\n', subj, processedSubsThisJob, totalSubsThisJob);
end

fprintf('Job %d completed processing %d subjects.\n', jobID, processedSubsThisJob);

%% Helper function to append to log file
function appendToLog(filename, message)
    fid = fopen(filename, 'a');
    if fid == -1
        warning('Could not open log file: %s', filename);
        return;
    end
    fprintf(fid, message);
    fclose(fid);
end