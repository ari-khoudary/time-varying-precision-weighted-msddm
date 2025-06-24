% compute slope of DV from pre-generated data structure

[allData(:).dvSlopes] = deal(zeros(100,6));

for i = 1:height(allData)
    % initialize variables
    noise1 = allData(i).noise1Frames;
    signal1 = noise1 + allData(i).signal1Frames;
    noise2 = noise1 + signal1 + allData(i).noise2Frames;
    signal2 = noise1 + signal1 + noise2 + allData(i).signal2Frames;
    dv = allData(i).decisionVariable;
    nFrames = allData.nFrames;
    RT = allData(i).RT;

    % compute slopes: (y2-y1) / (x2-x1)
    % first noise period
    allData(i).dvSlopes(:,1) = dv(noise1) ./ noise1;
    
    % first signal period
    allData(i).dvSlopes(:,2) = (dv(signal1) - dv(noise1)) ./ (signal1 - (noise1));

    % second noise period
    allData(i).dvSlopes(:,3) = (dv(noise2) - dv(signal1)) ./ (noise2 - (signal1));

    % second signal period
    allData(i).dvSlopes(:,4) = (dv(signal2) - dv(noise2)) ./ (signal2 - (noise2));

    % total slope change 
    allData(i).dvSlopes(:,5) = (dv(nFrames,:) - dv(1,:)) ./ nFrames;

    % total change w/ thresholded dv
    allData(i).dvSlopes(:,6) = (dv(RT) - dv(ones(size(RT)))) ./ nFrames;

end

clear noise1 noise2 signal1 signal2 dv RT

