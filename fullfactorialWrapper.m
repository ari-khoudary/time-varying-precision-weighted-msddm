% full factorial sim

cues = [0.5, 0.65, 0.8];
antCohs = [0.5, 0.65, 0.8];
coherences = [0.5, 0.65, 0.8];
congruence = [0, 1];

for i=1:length(cues)
    for j = 1:length(antCohs)
        for k = 1:length(coherences)
            for l = 1:length(congruence)
                doSampling(cues(i), antCohs(j), coherences(k), congruence(l));
            end
        end
    end
end
