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

# III. From posterior samples:
## A. Results from simulations: SYLVIA
### 1. Alphas (low prev), able to capture within- and among-patient corrs that are different
### 2. With fixed effects, recovery of all params
## B. Results from posterior samples: JOE
### 1. Benefits of using stan 
#### (solving some computational challenges (speed, efficiency, parallelizable), extendable to microbial data)
### 1. Alphas low
### 2. TBV effects on gam and phi
### 3. Correlations among and within patients
### 4. Fixed effects on phi and gamma
### 6. Limitations:
#### a. no assumptions about dynamics occurring between patient observations
#### b. not able to comment on the mechanisms leading to observed residual correlations
#### c. pairwise fixed effects should be interpretted with caution

# IV. Future extensions:
## A. Detection error
## B. blah
################################################
################################################