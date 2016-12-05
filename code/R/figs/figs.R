# Preliminary Figures 

################################################
################################################
# To Do:

# I. Model Motivation: Using a 2 species schematic
## A. How to quantify community dynamics through time? 
### 1. Data summaries and null models might obscure dynamics.
### 2. Species can be attracted to similar resources
### 3. Species can interact
### 4. These dynamics can be opposed at different scales
## B. Correlated random effects account for unknown shared resources or other unexplained dynamics
## C. Pairwise fixed effects allow for estimation of effects on colonization and extinction probabilities

# Required Schematics:
# 1. Ability to capture prevalence among sites
# 2. Potential dynamics:
## A. Habitat filtering, no species interactions (positive correlations, no fixed effects)
## B. Species interactions, no filtering (no correlations, only sig. fixed effects)
### Example: 2 species never occur in same hosts, because hosts are differnt, 
###          or they never occur due to interactions among hosts
### IMPORTANT: Either A or B could lead to these outcomes

## C. Temporal dynamics:
### 1. Mimic what's on the white board
### 2. Change the different questions as we change what's in the person boxes
### - Colonization:
#      - Positive effect, negative effect
### - Persistence:
#      - Positive effect, negative effect
#   3. Different patient types, different temporal dynamics (or same)
#   4. Patient-level and fixed effects could be different direction 

# II. From Raw Data: SYLVIA
## A. Rank-order prevalence (over time)
## B. Observed co-infect to single-infect ratios

# III. From postserior samples:
## A. Results from simulations: SYLVIA
### 1. Alphas (low prev), able to capture within- and among-patient corrs that are different
### 2. With fixed effects, recovery of all params
## B. Results from postserior samples: JOE
### 1. Benefits of using stan 
#### (solving some computational challenges (speed, efficiency, parallelizable), extendable to microbial data)
### 2. Alphas low
### 3. TBV effects on gam and phi
### 4. Correlations among and within patients
### 5. Fixed effects on phi and gamma
### 6. Limitations:
#### a. no assumptions about dynamics occurring between patient observations
#### b. not able to comment on the mechanisms leading to observed residual correlations
#### c. pairwise fixed effects should be interpretted with caution

# IV. Future extensions:
## A. Detection error
## B. blah
################################################
################################################
library(tidyverse)
library(scales)
library(ggsci)
library(rstan)
library(grid)
library(gridExtra)

# Do you want to save the figures?
save_it = 1
################################################
################################################

################################################
# Co-occurrence data from raw:
################################################

load(.) # from raw_data_analysis

df_cooccur = full_data_cooccur$results
#p_lt prob less than expected (negative)
#p_gt prob greater than expected (positive)

d = full_data_cooccur$species

df_cooccur = 
  df_cooccur %>%
  rowwise() %>%
  mutate(sig = ifelse(p_lt <= 0.05 | p_gt <= 0.05, 1, 0)) %>%
  mutate(sig = ifelse(sig == 1 & p_lt < 0.05, 2, sig)) 


plot_cooccur = 
  ggplot(df_cooccur, aes(x = sp1, y = sp2, fill = factor(sig))) +
  scale_fill_manual("Co-\noccurrence", values = c("gray", "blue", "red"),
                       labels = c("random", "positive", "negative"),
                       na.value = "white") +
  theme_classic() +
  scale_y_continuous(breaks = 1:d,
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     trans = "reverse") +
  scale_x_continuous(breaks = 1:d,
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     position = "top") +
  geom_tile() +
  labs(x="", y="") + 
  theme(axis.text.x.top = element_text(angle=30, vjust = 0.5),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

quartz(height=4, width=5, dpi=300)
plot_cooccur
if(save_it) ggsave("./figs/cooccur.pdf")

################################################
# Co-occurrence (obs v expect) data from raw:
################################################

plot_exp_cooccur = 
  ggplot(df_cooccur, aes(x = exp_cooccur, y = obs_cooccur, color = factor(sig))) +
  geom_point(shape = 19, size = 4) +
  scale_color_manual("Co-\noccurrence", values = c("gray", "blue", "red"),
                    labels = c("random", "positive", "negative"),
                    na.value = "gray") +
  scale_y_continuous(limits = c(0, 41), breaks = c(0, 20, 40)) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 15, 30)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_classic() +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black")) + 
  labs(x = "Expected Co-occurrences", y = "Observed Co-occurrences")

quartz(height=4, width=5, dpi=300)
plot_exp_cooccur
if(save_it) ggsave("./figs/exp_cooccur.pdf")
################################################
# Prevalence over time, raw data:
################################################

# dfm, from raw_data_analysis.R
dfm_last = filter(dfm, variable == "v10")

plot_prev =
  ggplot(dfm, aes(x = factor(variable), y = value, group = strain)) + 
  geom_line(aes(color = strain), size=1.25) + 
  scale_color_futurama(name = "Strain") + 
  scale_y_continuous(limits = c(0, 0.095), breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  scale_x_discrete(labels = 1:10) +
  annotate("text", x = dfm_last$variable, y = dfm_last$value, label = dfm_last$strain,
           size = 4) +
  theme_classic() +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, color = "black")) + 
  labs(x = "Visit", y = "Prevalence")

quartz(height=4, width=5.5, dpi = 300)
plot_prev
if(save_it) ggsave("./figs/prevalence.pdf")

################################################
# Simpson's paradox, from simulation:
################################################
# I manipulated the simulation script (sim_simpson.R) to produce:
# 1. Positive among host correlation
# 2. Negative among observation correlation

# Load the data I generated:
load("./code/R/figs/simpsons.RData")

# Store and plot it:
df_simpson <- data.frame(y_star, patient, visit)
colnames(df_simpson) = c("Strain1", "Strain2", "patient", "visit")

plot_simpson1 = 
  ggplot(df_simpson, aes(x=(Strain1), y=(Strain2))) +
  #ggplot(df_simpson, aes(x=pnorm(Strain1), y=pnorm(Strain2))) +
  geom_point(shape=19, aes(color=factor(patient))) +
  scale_color_futurama(name="Site", guide = F) +
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,0,3)) +
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,0,3)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,.5,1)) +
  #scale_x_continuous(limits = c(0,1), breaks = c(0,.5,1)) +
  theme_classic() +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold"))+
  labs(x = NULL, y = NULL)

quartz(height=5, width=5)
plot_simpson1
#if(save_it) ggsave("./figs/simpson1.pdf")

plot_simpson2 = 
  plot_simpson1 +
  facet_wrap(~patient, nrow = 2) +
  labs(x = NULL, y = NULL) +
  scale_color_futurama(name="Site") + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(1.2, "lines"))
  
quartz(height=5, width=10)
plot_simpson2

# Combine the figures:
plot_simpson =
  arrangeGrob(plot_simpson1, plot_simpson2,
             ncol = 3,
             layout_matrix = matrix(c(1,2,2), ncol = 3),
             left = textGrob("Occurrence Probability, Spp. 2", rot = 90,
                             gp=gpar(fontsize = 14, fontface = "bold")),
             bottom = textGrob("Occurrence Probability, Spp. 1", 
                               gp=gpar(fontsize = 14, fontface = "bold")))
quartz(height=5, width = 12)
grid.arrange(plot_simpson)
if(save_it) ggsave("./figs/simpson.pdf", plot_simpson, height=5, width=12, dpi=300)

################################################
# Results from postseriors:
################################################

# import the data:
load("./output/fit_full_HIM_10_strain.rda")

# extract postseriors
posts = extract(m_fit)

# dimensions (n_strain):
d = m_fit@par_dims$alphas
n_sample = dim(posts$alphas)[1]

############
############
# Alphas:
# Restructure for ggplot
alpha_post = posts$alphas
dim(alpha_post) = c(n_sample * d, 1)
alpha_post = data.frame(alpha_post)
colnames(alpha_post) = c("value")
alpha_post$strain = rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample)

plot_alphas =
  ggplot(alpha_post, aes(x=pnorm(value))) +
  geom_histogram(binwidth=0.0015, linetype = 0, fill="black") + 
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        strip.text = element_text(size = 12, face = "bold", color = "black")) + 
  facet_wrap(~strain, nrow=2, scales="free_x") +
  scale_y_continuous(limits=c(0,1800), breaks=c(0,900,1800)) +
  scale_x_continuous(limits=c(0, 0.1), breaks=seq(0,.08,.04)) +
  labs(x = expression(bold("Baseline occurrence probability, ") ~ bold(alpha)),
       y = "Frequency") 

quartz(height=5, width=10, dpi=300)
plot_alphas
if(save_it) ggsave("./figs/alphas.pdf")

############
############

# TBV effects on phi and gamma
tbv_phi = posts$betas_tbv_phi
tbv_gam = posts$betas_tbv_gam
dim(tbv_phi) = c(n_sample * d, 1)
dim(tbv_gam) = c(n_sample * d, 1)

tbv = data.frame(rbind(tbv_phi, tbv_gam))
colnames(tbv) = "value"
tbv$strain = rep( rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample), 2)
tbv$param = rep(c("phi", "gamma"), each = n_sample*d)

# which strains are significantly differnt than zero?
tbv_sigs = tbv %>%
  group_by(strain, param) %>%
  summarize(high = quantile(value, probs=c(0.975)),
            low = quantile(value, probs=c(0.025)),
            median = quantile(value, probs=c(0.5))) %>%
  mutate(sign_high = sign(high),
         sign_low = sign(low)) %>%
  rowwise() %>%
  mutate(sig = ifelse(sign_high*sign_low == 1, 1, 0)) %>%
  mutate(sig_pos = ifelse(sig == 1 & median > 0, 1, 0)) %>%
  filter(sig == 1) %>% 
  select(strain, param, sig_pos)

tbv$sig = 0
for(i in 1:nrow(tbv_sigs)){
  
  val = ifelse(tbv_sigs$sig_pos[i] == 0, 1, 2)
  tbv[tbv$strain == tbv_sigs$strain[i] & tbv$param == tbv_sigs$param[i], 4] = val
  
}

plot_tbv =
  ggplot(tbv, aes(x = value, fill=factor(sig))) +
  geom_histogram(binwidth=0.025, linetype=0) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        strip.text.y = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        strip.text.x = element_text(size = 12, face = "bold", color = "black")) + 
  guides(fill=F) + 
  facet_grid(param~strain, 
             labeller = labeller(param = label_parsed),
             scales = "free_x") +
  scale_y_continuous(limits=c(0,2000), breaks=c(0,1000,2000)) +
  scale_x_continuous(limits=c(-.6, .3), breaks = c(-.5, 0, 0.3)) +
  scale_fill_manual(values = c("black", "red", "blue")) +
  labs(x = expression(bold(beta[TBV])), y = "Frequency")
  
quartz(height=4, width=11, dpi=300)
plot_tbv
if(save_it) ggsave("./figs/tbv.pdf")

############
############

# Correlations within patients
corr_pat = posts$Rho_patient
dim(corr_pat) = c(n_sample*d*d,1)
colnames(corr_pat) = "value"
corr_pat = data.frame(corr_pat)
corr_pat$strain1 = rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample*d)
corr_pat$strain2 = rep(rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample), d)

# Summarize (median), and if not significant make NA
corr_pat_med = 
  corr_pat %>%
  group_by(strain1, strain2) %>%
  summarize(high = quantile(value, probs=c(0.975)),
            low = quantile(value, probs=c(0.025)),
            median = quantile(value, probs=c(0.5))) %>%
  mutate(sign_high = sign(high),
         sign_low = sign(low)) %>%
  rowwise() %>%
  mutate(sig = ifelse(sign_high*sign_low == 1, 1, 0)) %>%
  mutate(sig_pos = ifelse(sig == 1 & median > 0, 1, 0)) %>%
  select(strain1, strain2, median, sig, sig_pos) %>%
  mutate(med_fixed = ifelse(sig != 1 | median == 1.0, NA, median)) %>%
  mutate(med_fixed = round(med_fixed, 2))

# just upper tri
test_mat = matrix(1:100, nrow = 10, ncol = 10)
these = test_mat[lower.tri(test_mat)]

corr_pat_med = corr_pat_med[-these, ]
corr_pat_med$strain1 = as.numeric(as.factor((corr_pat_med$strain1)))
corr_pat_med$strain2 = as.numeric(as.factor((corr_pat_med$strain2)))

plot_corr_pat = 
  ggplot(corr_pat_med, aes(x = strain2, y = strain1, fill = med_fixed)) +
  geom_tile() +
  scale_fill_gradient2("Patient\nCorr", limits = c(-.5,.5), breaks = c(-.5,0,.5),
                      labels = c("-0.5", "  0", "  0.5"),
                      na.value = "white", low="red", mid = "white", high="blue") +
  theme_classic() +
  scale_y_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), 
                     trans="reverse") +
  scale_x_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     position = "top") +
  labs(x="", y="") + 
  theme(axis.text.x.top = element_text(angle=30, vjust = 0.5),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

quartz(height=4, width=5, dpi=300)
plot_corr_pat
if(save_it) ggsave("./figs/corr_pat.pdf")

############
############

# Correlations among observations
corr_obs = posts$Rho_visit
dim(corr_obs) = c(n_sample*d*d,1)
colnames(corr_obs) = "value"
corr_obs = data.frame(corr_obs)
corr_obs$strain1 = rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample*d)
corr_obs$strain2 = rep(rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample), d)

# Summarize (median), and if not significant make NA
corr_obs_med = 
  corr_obs %>%
  group_by(strain1, strain2) %>%
  summarize(high = quantile(value, probs=c(0.975)),
            low = quantile(value, probs=c(0.025)),
            median = quantile(value, probs=c(0.5))) %>%
  mutate(sign_high = sign(high),
         sign_low = sign(low)) %>%
  rowwise() %>%
  mutate(sig = ifelse(sign_high*sign_low == 1, 1, 0)) %>%
  mutate(sig_pos = ifelse(sig == 1 & median > 0, 1, 0)) %>%
  select(strain1, strain2, median, sig, sig_pos) %>%
  mutate(med_fixed = ifelse(sig != 1 | median == 1.0, NA, median)) %>%
  mutate(med_fixed = round(med_fixed, 2))

# just upper tri
test_mat = matrix(1:100, nrow = 10, ncol = 10)
these = test_mat[lower.tri(test_mat)]

corr_obs_med = corr_obs_med[-these, ]
corr_obs_med$strain1 = as.numeric(as.factor((corr_obs_med$strain1)))
corr_obs_med$strain2 = as.numeric(as.factor((corr_obs_med$strain2)))

plot_corr_obs = 
  ggplot(corr_obs_med, aes(x = strain2, y = strain1, fill = med_fixed)) +
  geom_tile() +
  scale_fill_gradient2("Obs\nCorr", limits = c(-.5,.5), breaks = c(-.5,0,.5),
                       labels = c("-0.5", "  0", "  0.5"),
                       na.value = "white", low="red", mid = "white", high="blue") +
  theme_classic() +
  scale_y_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), 
                     trans="reverse") +
  scale_x_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     position = "top") +
  labs(x="", y="") + 
  theme(axis.text.x.top = element_text(angle=30, vjust = 0.5),
        axis.text = element_text(size = 12, face = "bold", color = "black"))

quartz(height=4, width=5, dpi=300)
plot_corr_obs
if(save_it) ggsave("./figs/corr_obs.pdf")

############
############

# Fixed effects: Beta_phi
beta_phi = posts$betas_phi
dim(beta_phi) = c(n_sample*d*d,1)
colnames(beta_phi) = "value"
beta_phi = data.frame(beta_phi)
beta_phi$strain1 = rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample*d)
beta_phi$strain2 = rep(rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample), d)

# Summarize (median), and if not significant make NA
beta_phi_med = 
  beta_phi %>%
  group_by(strain1, strain2) %>%
  summarize(high = quantile(value, probs=c(0.975)),
            low = quantile(value, probs=c(0.025)),
            median = quantile(value, probs=c(0.5))) %>%
  mutate(sign_high = sign(high),
         sign_low = sign(low)) %>%
  rowwise() %>%
  mutate(sig = ifelse(sign_high*sign_low == 1, 1, 0)) %>%
  mutate(sig_pos = ifelse(sig == 1 & median > 0, 1, 0)) %>%
  select(strain1, strain2, median, sig, sig_pos) %>%
  mutate(med_fixed = ifelse(sig != 1 | median == 1.0, NA, median)) %>%
  mutate(med_fixed = round(med_fixed, 2))

beta_phi_med$strain1 = as.numeric(as.factor((beta_phi_med$strain1)))
beta_phi_med$strain2 = as.numeric(as.factor((beta_phi_med$strain2)))

plot_beta_phi = 
  ggplot(beta_phi_med, aes(x = strain1, y = strain2, fill = med_fixed)) +
  geom_tile() +
  scale_fill_gradient2(expression(beta[phi]), limits = c(-1,1), breaks = c(-1,0,1),
                       labels = c("-1", "  0", "  1"),
                       na.value = "white", low="red", mid = "white", high="blue") +
  theme_classic() +
  scale_y_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), 
                     trans="reverse") +
  scale_x_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     position = "top") +
  labs(x="", y="") + 
  theme(axis.text.x.top = element_text(angle=30, vjust = 0.5))

quartz(height=6, width=6, dpi=300)
plot_beta_phi
if(save_it) ggsave("./figs/beta_phi.pdf")

############
############

# Fixed effects: Beta_gam
beta_gam = posts$betas_gam
dim(beta_gam) = c(n_sample*d*d,1)
colnames(beta_gam) = "value"
beta_gam = data.frame(beta_gam)
beta_gam$strain1 = rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample*d)
beta_gam$strain2 = rep(rep(paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), each = n_sample), d)

# Summarize (median), and if not significant make NA
beta_gam_med = 
  beta_gam %>%
  group_by(strain1, strain2) %>%
  summarize(high = quantile(value, probs=c(0.975)),
            low = quantile(value, probs=c(0.025)),
            median = quantile(value, probs=c(0.5))) %>%
  mutate(sign_high = sign(high),
         sign_low = sign(low)) %>%
  rowwise() %>%
  mutate(sig = ifelse(sign_high*sign_low == 1, 1, 0)) %>%
  mutate(sig_pos = ifelse(sig == 1 & median > 0, 1, 0)) %>%
  select(strain1, strain2, median, sig, sig_pos) %>%
  mutate(med_fixed = ifelse(sig != 1 | median == 0.0, NA, median)) %>%
  mutate(med_fixed = round(med_fixed, 1))

beta_gam_med$strain1 = as.numeric(as.factor((beta_gam_med$strain1)))
beta_gam_med$strain2 = as.numeric(as.factor((beta_gam_med$strain2)))

plot_beta_gam = 
  ggplot(beta_gam_med, aes(x = strain1, y = strain2, fill = med_fixed)) +
  geom_tile() +
  scale_fill_gradient2(expression(beta[gamma]), limits = c(-1,1), breaks = c(-1,0,1),
                       labels = c("-1", "  0", "  1"),
                       na.value = "white", low="red", mid = "white", high="blue") +
  theme_classic() +
  scale_y_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)), 
                     trans="reverse") +
  scale_x_continuous(breaks = 1:d, 
                     labels = paste0("hpv",c(6,11,16,18,31,33,45,52,58,84)),
                     position = "top") +
  labs(x="", y="") + 
  theme(axis.text.x.top = element_text(angle=30, vjust = 0.5))

quartz(height=6, width=6, dpi=300)
plot_beta_gam
if(save_it) ggsave("./figs/beta_gam.pdf")

