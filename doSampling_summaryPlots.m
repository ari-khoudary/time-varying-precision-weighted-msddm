%% doSampling_summaryPlots

%% load in .mat files
dataDir = outDir;
files = (dir([pwd filesep dataDir filesep '*.mat']));
allData = load([outDir filesep files(1).name]).data;
allData = repmat(allData, length(files), 1);
for i=2:length(files)
    allData(i) = load([outDir filesep files(i).name]).data;
end

%% add fields to be populated

for i = 1:length(allData)
    congruent = logical(allData(i).congruent);
    allData(i).congruentRawAccuracy = mean(allData(i).choices(congruent, 1));
    allData(i).congruentForcedAccuracy = mean(allData(i).choices(congruent, 2));
    allData(i).incongruentRawAccuracy = mean(allData(i).choices(~congruent, 1));
    allData(i).incongruentForcedAccuracy = mean(allData(i).choices(~congruent,2));
    allData(i).congruentRT = mean(allData(i).RT(congruent));
    allData(i).incongruentRT = mean(allData(i).RT(~congruent));
end






