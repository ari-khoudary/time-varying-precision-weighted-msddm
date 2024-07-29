# trial plots
library(tidyverse)
blue = 'dodgerblue1'
orange = 'darkorange'
yellow = 'sienna4'

late = read.csv('2022_ccn/poster/0.80cue_0.50coh_trial6.csv')
early = read.csv('2022_ccn/poster/0.80cue_0.50coh_trial8.csv')
threshold=24

## late effect
late <- late %>%
  mutate(time = -10:99) %>% 
  rename(memoryEvidence = mem,
         visualEvidence = viz)
late_idx = min(which(late$dv > threshold))

# trace
late[11:late_idx,] %>%
  ggplot(aes(x=time)) +
  geom_hline(yintercept = c(0,-25,25), linetype=c('dotted', 'solid','solid'), linewidth=c(0.5, 0.25, 0.25)) + 
  geom_vline(xintercept = 0, linetype='dashed', linewidth=0.25) + 
  geom_line(aes(y=memoryEvidence), color=blue, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=visualEvidence), color = orange, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=dv), color = yellow, linewidth=1.5) +
  #geom_point(aes(y=visualEvidence[11], x=0),color=orange, size=1.5) +
  #geom_point(aes(y=dv[11], x=0),color=yellow, size=1.5) +
  ylim(-30, 100) + 
  xlim(-10, 90) + 
  theme_void() 
ggsave('lateEffect.png', width=4, height=3)

        
### early effect
early <- early %>%
  mutate(time = -10:99) %>%
  rename(memoryEvidence = mem,
         visualEvidence = viz)
early_idx = min(which(early$dv > threshold))

# trace
early[1:early_idx,] %>%
  ggplot(aes(x=time)) +
  geom_hline(yintercept = c(0,-25,25), linetype=c('dotted', 'solid','solid'), linewidth=c(0.5, 0.25, 0.25)) + 
  geom_vline(xintercept = 0, linetype='dashed', linewidth=0.25) + 
  geom_line(aes(y=memoryEvidence), color=blue, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=visualEvidence), color = orange, linewidth=1, alpha=0.7) + 
  geom_line(aes(y=dv), color = yellow, linewidth=1.5) +
  #geom_point(aes(y=dv[11], x=0), color=yellow, size=1.5) +
  #geom_point(aes(y=visualEvidence[11], x=0), color=orange, size=1.5) +
  ylim(-30, 100) + 
  xlim(-10, 90) + 
  theme_void()
ggsave('earlyEffect.png', width=4, height=3)


######### drift value plots #########
# late
ggplot(late, aes(x=time)) +
  geom_vline(xintercept = 0, linetype='dashed') + 
  geom_line(aes(y=memoryDrift), color=blue, size=2) + 
  geom_line(aes(y=visualDrift), color = orange, size=2) + 
  xlim(-20, 90) + 
  theme_classic() +
  theme(text = element_text(size=14, family='Avenir')) + 
  labs(y='weight (a.u.)', x='time (a.u.)')
ggsave('lateDrift.png', width=6, height=2)

# early
ggplot(early, aes(x=time)) +
  geom_vline(xintercept = 0, linetype='dashed') + 
  geom_line(aes(y=memoryDrift), color=blue, size=2) + 
  geom_line(aes(y=visualDrift), color = orange, size=2) + 
  scale_x_continuous(breaks = c(-15, 0, 30, 60, 90), 
                     labels = c(-10, 0, 100, 200, 300),
                     limits = c(-20, 100)) + 
  theme_classic() +
  theme(text = element_text(size=24, family='Avenir')) + 
  labs(y='drift rate')
ggsave('earlyDrift.png', width=6, height=4)
