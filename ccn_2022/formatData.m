%% make csvs for plotting

files = dir('results/*.mat');
for i=1:length(files)
    % load up each .mat
    infile = files(i).name; 
    simID=char(extractBefore(infile, '_100subs'));
    load(infile);

    % get integrated evidence accumulation matrix for each subject
    fullEvidence_long = permute(fullEvidence, [3 2 1]);
    
    for subj=1:size(fullEvidence_long,3)
        % make each subject a csv, add their subID
        evidenceOutfile = fullEvidence_long(:, :, subj);
        evidenceOutfile(:, end+1) = subj;
    
        % make rows of behavior = trials, add subID
        behaviorOutfile = choices(subj, :)';
        behaviorOutfile(:,end+1) = RT(subj, :)';
        behaviorOutfile(:, end+1) = subj;
    
        % write files
        evidenceFilename = sprintf('results/evidenceMatrices/%s_s%i.csv', simID, subj);
        behaviorFilename = sprintf('results/behavior/%s_s%i.csv', simID, subj);
        csvwrite(evidenceFilename, evidenceOutfile);
        csvwrite(behaviorFilename, behaviorOutfile);
    end
end