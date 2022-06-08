# Precision-weighted evidence integration predicts time-varying influence of memory on perceptual decisions

*Abstract*: How do we use past experiences to make sense of the present? We developed an experimental task and computational model to investigate whether human observers Bayes-optimally integrate evidence from memory and vision during cued perceptual decisions. Drawing on theory developed in the multisensory integration literature, we model the decision process as precision-weighted integration of samples from memory and vision, with the rate of accumulation defined by the relative precision of each evidence stream. The model captures two qualitatively distinct empirical findings from different task structures, suggesting this framework has the potential to identify a process fundamental to decisions relying on multiple evidence sources. The experimental task will measure how human observers actually perform this integration, which will allow us to better characterize the conditions under which humans deviate from optimal choice behavior.

# About this repo
- `doSampling.m` generates samples, timecourses, decisions, and reaction times for `nSub` completing `nTrial`. It stores the output of each simulation in `results`. 
- `plotTraces.m` plots timecourses for memory, visual, and combined sampling for a given results file
- `plotDrift.m` plots the time-varying drift rates of memory and visual evidence for a given results file
- `plotPrecision.m` plots the time-varying precision estimates of memory and visual evidence for a given results file.

# Questions? Comments?
Email makhouda [at] uci.ecu
