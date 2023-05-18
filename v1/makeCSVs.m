dir = 'bayes_results/';
cues = [0.50, 0.65, 0.8];
antCohs = [0.50, 0.65, 0.8];
coherences = [0.50, 0.65, 0.8];
congruence = [0, 1];
v = [10, 1000];

for cue=1:length(cues)
    for j = 1:length(antCohs)
        for k = 1:length(coherences)
            for l = 1:length(congruence)
                for m = 1:length(v)
                % load in the file
                infile=sprintf('%s%.2fcue_%.2fantCoh_%.2fcoh_%icong_25thresh_4memThin_altBeta_%iv.mat', dir, cues(cue), antCohs(j), coherences(k), congruence(l), v(m));
                load(infile);
                accuracy = choices';
                RT = RT';
                cueLevel = repmat(cueLevel, [nTrial 1]);
                anticipatedCoherence = repmat(anticipatedCoherence, [nTrial 1]);
                coherenceLevel = repmat(coherenceLevel, [nTrial 1]);
                congruent = repmat(congruent, [nTrial 1]);
                v = repmat(v, [nTrial 1]);
                behav_csv = table(accuracy, RT, cueLevel, anticipatedCoherence, coherenceLevel, congruent, v);
                filename = [dir 'behav_csv/' char(extractBetween(outfile, '/', '.mat')) '.csv'];
                writetable(behav_csv, filename);
                end
            end
        end
    end
end