# plot alpha beta timecourses for poster traces
library(tidyverse)
library(ggdist)
library(distributional)
blue = 'dodgerblue1'
orange = 'darkorange1'
yellow = '#ffd103'

# load
bornstein_df = read.csv('bornsteinEffect.csv')
hanks_df = read.csv('hanksEffect.csv')

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
