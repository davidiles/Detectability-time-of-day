
# Simulate an effect of time-since-sunrise (TSS) that differs between each of two sexes
#  - then evaluate the ability of various approaches to correct for those effects and estimate true (or even just relative) density



male_proportion <- 0.5  # Proportion of males in population
n_surveys <- 1000

male_density <- 1 # birds per ha

survey_duration <- 3  # minutes
time_bins <- c(0,1,2)

# --------------------------------
# At each survey location, assume that 
# --------------------------------