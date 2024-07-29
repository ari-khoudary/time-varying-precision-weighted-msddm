% make csvs

%% load in data
b = load('bornstein_8thin.mat');
h = load('hanks_8thin.mat');
b = b.data;
h = h.data;

trial_b = 7;
trial_h = 3;

%% make tables

b_result = table(b.counters(:,1, trial_b), b.counters(:,2,trial_b), b.counters(:,3,trial_b), b.counters(:,4,trial_b), ...
    b.memoryAccumulator(:, trial_b), b.visionAccumulator(:,trial_b), b.decisionVariable(:,trial_b), [1:b.nFrames]', ...
    repelem(b.trialDelays(trial_b), b.nFrames)');
b_result = renamevars(b_result, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9"], ...
    ["memAlpha", "memBeta", "vizAlpha", "vizBeta", "memory", "vision", "DV", "timestep", "trialDelay"]);

h_result = table(h.counters(:,1, trial_h), h.counters(:,2,trial_h), h.counters(:,3,trial_h), h.counters(:,4,trial_h), ...
    h.memoryAccumulator(:, trial_h), h.visionAccumulator(:,trial_h), h.decisionVariable(:,trial_h), [1:h.nFrames]');
h_result = renamevars(h_result, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8"], ...
    ["memAlpha", "memBeta", "vizAlpha", "vizBeta", "memory", "vision", "DV", "timestep"]);

%% write csvs

writetable(b_result, 'bornsteinEffect.csv')
writetable(h_results, 'hanksEffect.csv')
