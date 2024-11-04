# Social network dynamics in response to parasitism (guppy-gyro model) #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 1 - Data wrangling



#### Setup ####

## Get the environment ready for loading the data

# Clear environment
rm(list=ls())

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Social networks in guppies (Scientific Reports)")

# Load the necessary libraries
library(igraph); library(ggplot2); library(plyr); library(dplyr)
library(stringr); library(purrr); library(gridExtra); library(reshape2)
library(lme4); library(lmerTest); library(MuMIn); library(DHARMa)
library(glmmTMB); library(car); library(ggraph); library(ggnetwork)

#### Data input ####

## Read in the data

# Read in social network data in an edgelist format
social_edges_long <- read.csv("Data/social_nets.csv")
social_edges <- subset(social_edges_long, shoal != "L")

# Infection data for shoals 
infections <- read.csv("Data/ind_infect.csv")

# Meta-data for infection start points
targets <- read.csv("Data/ind_target.csv")


#### Data manipulation and wrangling #### 

## Split into the different types of interaction

# Create the sum of interactions for each pair 
edges <- aggregate(. ~ shoal * time * lower * upper, 
                   data = dplyr::select(social_edges, -duration), sum)

# Parasite intensity for the lower nodes (i.e., Fish.ID)
lower_infections <- left_join(edges, infections, by = c("shoal" = "shoal", 
                                                        "time" = "time",
                                                        "lower" = "fish_id"))

# Parasite intensity for the lower nodes (i.e., Fish.ID)
upper_infections <- left_join(edges, infections, 
                              by = c("shoal" = "shoal", 
                                     "time" = "time",
                                     "upper" = "fish_id"))

# Add to edges dataframe
edges$lower_infection <- lower_infections$parasite_intensity
edges$upper_infection <- upper_infections$parasite_intensity

## Separate the data out into lists for each shoal (replicate)

# Split the direct interaction data out into shoals (replicates)
edges_list <- split(edges, list(edges$shoal, edges$time))

## Create networks from the edgelists 

# Direction interactions
nets_list <- lapply(edges_list, function(x){graph_from_data_frame(as.matrix(dplyr::select(x, lower, upper, weight)), directed = T)})


## Create averaged networks before and after infection for one shoal (I) for the 
## final figure 

# Subset shoal I
edges_I <- subset(edges, shoal == "I")

# Add before and after categories (before and after infections)
edges_I$infection <-  cut(as.numeric(edges_I$time), c(0,5,10))
levels(edges_I$infection) <- c("Before", "After")

# Create unique edgelists for the networks
edges_ba <- aggregate(weight ~ upper*lower + infection, data = edges_I, FUN = "sum")

# Remove self loop
edges_ba <- edges_ba[-1,]

# Create the networks for before and after 
edges_I_ba <- split(edges_ba, edges_ba$infection)

# Direction interactions
nets_I_ba <- lapply(edges_I_ba, function(x){graph_from_data_frame(as.matrix(dplyr::select(x, lower, upper, weight)), directed = T)})

# Infections
infections_I <- subset(infections, shoal == "I")

# Add before and after categories (before and after infections)
infections_I$infection <-  cut(as.numeric(infections_I$time), c(0,5,10))
levels(infections_I$infection) <- c("Before", "After")

# Create category for infected and uninfected
infections_I_ba <- aggregate(parasite_intensity~infection + fish_id, data = infections_I, FUN = "mean")
infections_I_ba$colour <-  cut(as.numeric(infections_I_ba$parasite_intensity), c(-Inf,0,Inf))
levels(infections_I_ba$colour) <- c("darkgrey", "darkred")

# Order the fish IDs to match the order in the networks 
infections_I_b <- subset(infections_I_ba, infection == "Before")
infections_I_b$fish_id <- factor(infections_I_b$fish_id, levels = V(nets_I_ba[[1]])$name)
para_I_b <- infections_I_b[order(infections_I_b$fish_id), ]

# Order the fish IDs to match the order in the networks 
infections_I_a <- subset(infections_I_ba, infection == "After")
infections_I_a$fish_id <- factor(infections_I_a$fish_id, levels = V(nets_I_ba[[2]])$name)
para_I_a <- infections_I_a[order(infections_I_a$fish_id), ]

# Add the colour data for the nodes
nets_I_b <- nets_I_ba[[1]]
V(nets_I_b)$infection <- para_I_b[, "colour"]

# Add the colour data for the nodes
nets_I_a <- nets_I_ba[[2]]
V(nets_I_a)$infection <- para_I_a[, "colour"]

# Remove objects that are not needed
rm(lower_infections, upper_infections, social_edges, social_edges_long)
