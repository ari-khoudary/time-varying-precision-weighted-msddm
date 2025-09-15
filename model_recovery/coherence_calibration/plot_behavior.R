library(tidyverse)
library(patchwork)

df <- read.csv('results/results_2025-08-15_18-03-56.csv')

summary_df <- df %>% group_by(coherence, threshold, subID) %>%
  mutate(rawChoiceAccuracy = ifelse(noChoice==1 & rawChoiceAccuracy==0, NA, rawChoiceAccuracy)) %>%
  summarise(rawChoice_accuracy = mean(rawChoiceAccuracy, na.rm=T),
            forcedChoice_accuracy = mean(forcedChoiceAccuracy),
            rawChoice_rt = mean(rawChoiceRT, na.rm=T))

write.csv(summary_df, 'results/results_2025-08-15_summary_naAcc.csv', row.names = F)

## plot accuracy: raw choice
summary_df %>%
  mutate(threshold = factor(threshold)) %>%
  ggplot(aes(x=coherence, y=rawChoice_accuracy, color=threshold, group=threshold)) +
  theme_bw() + geom_hline(yintercept = 0.5, color='gray30') +
  stat_summary(fun = 'mean', geom='line') +
  stat_summary(fun.data = 'mean_se', geom='pointrange') +
  geom_hline(yintercept = 0.7, linetype = 'dashed') + 
  facet_wrap(~ threshold) +
  labs(title = 'raw choice accuracy', subtitle = '75,000 observations per point') | 

## plot accuracy: forced choice
summary_df %>%
  mutate(threshold = factor(threshold)) %>%
  ggplot(aes(x=coherence, y=forcedChoice_accuracy, color=threshold, group=threshold)) +
  theme_bw() + geom_hline(yintercept = 0.5, color='gray30') +
  stat_summary(fun = 'mean', geom='line') +
  stat_summary(fun.data = 'mean_se', geom='pointrange') +
  geom_hline(yintercept = 0.7, linetype = 'dashed') + 
  facet_wrap(~ threshold) +
  labs(title = 'forced choice accuracy', subtitle = '75,000 observations per point') 
  

## plot RT: raw choice
summary_df %>%
  mutate(threshold = factor(threshold)) %>%
  ggplot(aes(x=coherence, y=rawChoice_rt, color=threshold, group=threshold)) +
  theme_bw() + #geom_hline(yintercept = 0.5, color='gray30') +
  stat_summary(fun = 'mean', geom='line') +
  stat_summary(fun.data = 'mean_se', geom='pointrange') +
  # geom_hline(yintercept = 0.7, linetype = 'dashed') + 
  facet_wrap(~ threshold) +
  labs(title = 'raw choice RT', subtitle = '75,000 observations per point') 

