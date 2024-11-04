# Social network dynamics in response to parasitism (guppy-gyro model) #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 4 - Statistical analysis


#### Setup ####

## Get the environment ready for loading the data

# Clear environment
rm(list=ls())

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Social networks in guppies (Scientific Reports)")

# Load the previous scripts 
source("Code/3_Network-metrics.R")


#### Node metrics #### 

## Prepare the data 

# Scale parasite intensity (as there is a non-normal distribution of values 
# affecting the model fit)
inf_node_metrics$para_stand <- scale(inf_node_metrics$parasite_intensity)

# Trim the data to records with parasite intensity measured
inf_node_metrics_model <- subset(inf_node_metrics, parasite_intensity != "NA")

## Weighted in degree (poisson GLMM)

# Model structure
model1a <- glmer.nb(windegree ~ para_stand + 
                      treatment * infection + 
                      (1|shoal) + 
                      (1|fish_id),
                 data = inf_node_metrics_model)

# Model validation
simulationOutput1a <- simulateResiduals(fittedModel = model1a, plot = T)
testOutliers(simulationOutput = simulationOutput1a, type = "bootstrap")
testZeroInflation(simulationOutput1a)

# Null model
model1a_null <- glmer.nb(windegree ~ 1 + (1|shoal) + (1|fish_id), 
                                    data = inf_node_metrics_model)

# Model results
summary(model1a)
Anova(model1a)
anova(model1a_null, model1a)
r.squaredGLMM(model1a)


## Weighted out degree (Zero-inflated poisson GLMM)

# Model structure
model1b <- glmmTMB(woutdegree ~ para_stand + 
                      treatment * infection +
                      (1|shoal) + 
                      (1|fish_id),
                   ziformula = ~1, 
                 family = genpois(link = "log"),
                 data = inf_node_metrics_model)

# Model validation
simulationOutput1b <- simulateResiduals(fittedModel = model1b, plot = T)
testOutliers(simulationOutput = simulationOutput1b, type = "bootstrap")

# Null model
model1b_null <- glmmTMB(woutdegree ~ 1 + (1|shoal) + (1|fish_id),
                        #ziformula = ~1, 
                        family = genpois(link = "log"),
                        data = inf_node_metrics_model)

# Model results
summary(model1b)
car::Anova(model1b)
anova(model1b_null, model1b)
r.squaredGLMM(model1b)

# This is not the best model but I am struggling to troubleshoot - 


## Betweenness (Zero-inflated poisson GLMM)

# Model structure
model1c <- glmmTMB(round(betweenness) ~ para_stand + 
                   treatment * infection +
                   (1|shoal) + 
                   (1|fish_id),
                   ziformula = ~1, 
                   family = genpois(link = "log"),
                 data = inf_node_metrics_model)

# Model validation
simulationOutput1c <- simulateResiduals(fittedModel = model1c, plot = T)
testOutliers(simulationOutput = simulationOutput1c, type = "bootstrap")

# Null model
model1c_null <- glmmTMB(round(betweenness) ~ 1 + (1|shoal) + (1|fish_id),
                        #ziformula = ~1, 
                        family = genpois(link = "log"),
                        data = inf_node_metrics_model)

# Model results
summary(model1c)
Anova(model1c)
anova(model1c_null, model1c)
r.squaredGLMM(model1c)


#### Motifs #### 

## Data manipulation

# Split the data from the long dataframe 
motif_asymmetric_n <- subset(inf_motifs_long, variable == "inf_n")
motif_asymmetric_q <- subset(inf_motifs_long, variable == "inf_q")
motif_outstar_n <- subset(inf_motifs_long, variable == "os_n")
motif_outstar_q <- subset(inf_motifs_long, variable == "os_q")

## Asymmetric (count) GLMM 

# Create the relative motif frequency variable
motif_asymmetric_n$motif_s <- motif_asymmetric_n$value/motif_asymmetric_n$infected_fish_no

# Model structure
model2a <- lmer(motif_s ~ parasite_intensity +
                  treatment + 
                  #time + 
                  (1|shoal),
                 data = motif_asymmetric_n[-67,])

# Model validation
simulationOutput2a <- simulateResiduals(fittedModel = model2a, plot = T)
testOutliers(simulationOutput = simulationOutput2a, type = "bootstrap")
testZeroInflation(simulationOutput2a)

# Null model
model2a_null <- lmer(motif_s ~ 1 + (1|shoal),
                      data = motif_asymmetric_n[-67,])

# Model results
summary(model2a)
Anova(model2a)
anova(model2a_null, model2a)
r.squaredGLMM(model2a)


## Asymmetric (weighted count) GLMM 

# Create the relative motif frequency variable
motif_asymmetric_q$motif_s <- motif_asymmetric_q$value/motif_asymmetric_q$infected_fish_no

# Model structure
model2b <- lmer(motif_s ~ parasite_intensity +
                   treatment + 
                   #time + 
                   (1|shoal),
                 data = motif_asymmetric_q[-67,])

# Model validation
simulationOutput2b <- simulateResiduals(fittedModel = model2b, plot = T)
testOutliers(simulationOutput = simulationOutput2b, type = "bootstrap")
testZeroInflation(simulationOutput2b)

# Null model
model2b_null <- lmer(motif_s ~ 1 + (1|shoal),
                      data = motif_asymmetric_q[-67,])

# Model results
summary(model2b)
Anova(model2b)
anova(model2b_null, model2b)
r.squaredGLMM(model2b)


## Out-star (count) GLMM 

# Create the relative motif frequency variable
motif_outstar_n$motif_s <- motif_outstar_n$value/motif_outstar_n$infected_fish_no

# Model structure
model2c <- lmer(motif_s ~ parasite_intensity +
                   treatment + 
                   #time + 
                   (1|shoal),
                 data = motif_outstar_n[-67,])

# Model validation
simulationOutput2c <- simulateResiduals(fittedModel = model2c, plot = T)
testOutliers(simulationOutput = simulationOutput2c, type = "bootstrap")
testZeroInflation(simulationOutput2c)

# Null model
model2c_null <- lmer(motif_s ~ 1 + (1|shoal),
                      data = motif_outstar_n[-67,])

# Model results
summary(model2c)
Anova(model2c)
anova(model2c_null, model2c)
r.squaredGLMM(model2c)


## Out-star (weighted count) GLMM 

# Create the relative motif frequency variable
motif_outstar_q$motif_s <- motif_outstar_q$value/motif_outstar_q$infected_fish_no
motif_outstar_q$para_stand <- scale(motif_outstar_q$parasite_intensity)

# Model structure
model2d <- glmer(value ~ para_stand + 
                   infected_fish_no + 
                   #treatment + 
                   (1|shoal),
                 family = "poisson" (link = "log"),
                 data = motif_outstar_q)

# Model validation
simulationOutput2d <- simulateResiduals(fittedModel = model2d, plot = T)
testOutliers(simulationOutput = simulationOutput2d, type = "bootstrap")
testZeroInflation(simulationOutput2d)

# Null model
model2d_null <- glmer(value ~ 1 + (1|shoal), 
                      family = "poisson" (link = "log"),
                      data = motif_outstar_q)

# Model results
summary(model2d)
Anova(model2d)
anova(model2d_null, model2d)
r.squaredGLMM(model2d)


#### Network metrics #### 

## Connectance (edge density)

# Model structure
model3a <- glmer.nb(C ~ parasite_intensity +
                   treatment * infection +
                   (1|shoal),
                 #family = "binomial", 
                 data = inf_net_metrics)

# Model validation
simulationOutput3a <- simulateResiduals(fittedModel = model3a, plot = T)
testOutliers(simulationOutput = simulationOutput3a, type = "bootstrap")
testZeroInflation(simulationOutput3a)

# Null model
model3a_null <- glmer.nb(C ~ 1 + (1|shoal), 
                      data = inf_net_metrics)

# Model results
summary(model3a)
Anova(model3a)
anova(model3a_null, model3a)
r.squaredGLMM(model3a)


## Reciprocity 

# Model structure
model3b <- glmer.nb(reciprocity ~ parasite_intensity +
                  treatment * infection +
                  (1|shoal),
                data = inf_net_metrics)

# Model validation
simulationOutput3b <- simulateResiduals(fittedModel = model3b, plot = T)
testOutliers(simulationOutput = simulationOutput3b, type = "bootstrap")
testZeroInflation(simulationOutput3b)

# Null model
model3b_null <- lmer(C ~ 1 + (1|shoal), 
                     data = inf_net_metrics)

# Model results
summary(model3b)
Anova(model3b)
anova(model3b_null, model3b)
r.squaredGLMM(model3b)

