% simulate second-order regular Beta model on an ecologically valid scale
order = [2];
cues = [0.5, 0.65, 0.8];
antCohs = [0.5, 0.65, 0.8];
coherences = [0.5, 0.65, 0.8];
congruence = [0, 1];
beta = [0];
thresholds = [10, 25, 40];
memThin = [2, 4, 8];

tic
for i=1:length(order)
    for j=1:length(cues)
        for k = 1:length(antCohs)
            for l = 1:length(coherences)
                for m = 1:length(congruence)
                    for n = 1:length(beta)
                        for  o = 1:length(thresholds)
                            for  p = 1:length(memThin)
                                % don't need to pass "save" because i
                                % adjusted doSampling to auto-save within
                                % the subject loop
                                    doSampling(order(i), cues(j), antCohs(k), coherences(l), congruence(m), beta(n), thresholds(o), memThin(p), 0, 0, 0);
                            end
                        end
                     end
                end
            end
        end
    end
end
toc
