# plot alpha beta timecourses for poster traces
library(tidyverse)
library(ggdist)
library(distributional)
blue = 'dodgerblue1'
orange = 'darkorange1'
yellow = 'goldenrod'

# load
bornstein_df = read.csv('bornsteinEffect.csv')
hanks_df = read.csv('hanksEffect.csv')

# plot late effect
threshold = 10
late_idx = min(which(hanks_df$DV > threshold))

hanks_df[1:late_idx,] %>%
  ggplot(aes(x=timestep)) +
  geom_hline(yintercept = c(0,-threshold,threshold), linetype=c('dotted', 'solid','solid'), linewidth=c(0.5, 0.5, 0.5)) + 
  geom_vline(xintercept = 0, linetype='solid', linewidth=0.25) + 
  geom_line(aes(y=memory), color=blue, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=vision), color = orange, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=DV), color = yellow, linewidth=1.5) +
  labs(x = 'time (a.u.)', y='evidence') +
  theme_void()
ggsave('lateEffect.png', width=4, height=3, dpi='retina')

# plot early effect
threshold = 15
vizOnset = bornstein_df$trialDelay[1]
early_idx = min(which(bornstein_df$DV > threshold))
bornstein_df[1:early_idx,] %>%
  mutate(DV = ifelse(timestep<vizOnset, NA, DV),
         vision = ifelse(timestep<vizOnset, NA, vision)) %>%
  ggplot(aes(x=timestep)) +
  geom_hline(yintercept = c(0,-threshold,threshold), linetype=c('dotted', 'solid','solid'), linewidth=c(0.5, 0.5, 0.5)) + 
  geom_vline(xintercept = c(0, bornstein_df$trialDelay[1]), linetype=c('dashed', 'solid'), linewidth=c(0.25, 0.25)) + 
  geom_line(aes(y=memory), color=blue, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=vision), color = orange, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=DV), color = yellow, linewidth=1) +
  geom_point(aes(x=vizOnset, y=DV[vizOnset]), fill=yellow, color='gray20',size=2.5, shape=21) +
  geom_point(aes(x=vizOnset, y=vision[vizOnset]), fill=orange, color='gray20',size=2.5, shape=21) +
  labs(x = 'time (a.u.)', y='evidence') +
  theme_void()
ggsave('earlyEffect.png', width=4, height=3, dpi='retina')

# plot precision estimates for late effect
hanks_df[1:late_idx,] %>%
  ggplot(aes(x=timestep)) +
  geom_vline(xintercept = 0, linewidth=0.5, linetype='solid') + 
  geom_hline(yintercept = 1.5) +
  geom_line(aes(y=memPrecision), color=blue) +
  geom_line(aes(y=vizPrecision), color=orange) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = 'dynamic precision estimate (a.u.)', x = 'timestep (a.u.)')
ggsave('lateEffect_precisions.png', width=4, height=3, dpi='retina')

# plot precision estimates for early effect
bornstein_df[1:early_idx,] %>%
  mutate(vizPrecision=ifelse(timestep<vizOnset, NA, vizPrecision)) %>%
  ggplot(aes(x=timestep)) +
  geom_vline(xintercept = c(0, vizOnset), linewidth=0.5, linetype=c('dashed','solid')) + 
  geom_hline(yintercept = 1.5) +
  geom_line(aes(y=memPrecision), color=blue) +
  geom_line(aes(y=vizPrecision), color=orange) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = 'dynamic precision estimate (a.u.)', x = 'timestep (a.u.)')
ggsave('Effect_precisions.png', width=4, height=3, dpi='retina')

# plot
nDists = 20
subrows = round(seq(1, max(hanks_df$timestep), length.out=nDists))
hanks_df[subrows,] %>%
  ggplot(aes(x=timestep, ydist=dist_beta(memAlpha, memBeta))) +
  scale_x_continuous(breaks = subrows) +
  stat_slabinterval(color=blue, fill=blue, show_interval = FALSE) +
  stat_slabinterval(aes(ydist=dist_beta(vizAlpha, vizBeta)), color=orange, fill=orange, show_interval = FALSE) + 
  theme_classic() +
  labs(y = 'g_s(t): belief about value of generating parameter', title = 'late effect', x = 'timestep (~ms)')
ggsave('hanks_effect_betaPDFs.png', width=6, height=4, dpi='retina')


subrows = round(seq(1, max(bornstein_df$timestep), length.out=nDists))
bornstein_df[subrows,] %>%
  ggplot(aes(x=timestep, ydist=dist_beta(memAlpha, memBeta))) +
  scale_x_continuous(breaks = subrows) +
  stat_slabinterval(color=blue, fill=blue, show_interval = FALSE) +
  stat_slabinterval(aes(ydist=dist_beta(vizAlpha, vizBeta)), color=orange, fill=orange, show_interval = FALSE) + 
  theme_classic() +
  labs(y = 'g_s(t): belief about value of generating parameter', title = 'late effect', x = 'timestep (~ms)')
ggsave('bornstein_effect_betaPDFs.png', width=6, height=4, dpi='retina')
