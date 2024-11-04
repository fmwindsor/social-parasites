# Social network dynamics in response to parasitism (guppy-gyro model) #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 3 - Network metrics


#### Setup ####

## Get the environment ready for loading the data

# Clear environment
rm(list=ls())

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Social networks in guppies (Scientific Reports)")

# Load the previous scripts 
source("Code/1_Data-wrangling.R")

# Load the motif detector function
source("Code/Functions/Motif_detector.R")


#### Node metrics ####

## Create a dataframe to store the network metrics 

# Calculate weighted in degree
windegree <- unlist(lapply(nets_list, 
                           function(x){igraph::strength(x, mode = "in")}))

# Calculate weighted out degree
woutdegree <- unlist(lapply(nets_list, 
                            function(x){igraph::strength(x, mode = "out")}))

# Calculate betweenness for the inverse of interaction weights
betweenness <- unlist(lapply(nets_list, function(x){igraph::betweenness(x, weights = as.numeric(E(x)$weight), directed = T)}))

# Calculate closeness for the inverse of interaction weights
closeness <- unlist(lapply(nets_list, function(x){igraph::closeness(x, weights = as.numeric(E(x)$weight), normalized = T, mode = "total")}))

# Data frame for results
node_metrics <- data.frame(windegree, woutdegree, betweenness, closeness)

# Extract information from names of rows 
node_metrics$shoal <- unlist(map(str_split(rownames(node_metrics), "\\."), 1))
node_metrics$time <- unlist(map(str_split(rownames(node_metrics), "\\."), 2))
node_metrics$fish_id <- unlist(map(str_split(rownames(node_metrics), "\\."), 3))

# Add before and after categories (before and after infections)
node_metrics$infection <-  cut(as.numeric(node_metrics$time), c(0,5,10))
levels(node_metrics$infection) <- c("Before", "After")

# Add the treatments applied (infecting the most and least connected and a control)
node_metrics$treatment <-  revalue(node_metrics$shoal, 
                                         c("Aa" = "Least", "G" = "Least", 
                                           "H" = "Control", "I" = "Most",
                                           "J" = "Least", "K" = "Control", 
                                           "M" = "Most", "N" = "Least", 
                                           "O" = "Control", "P" = "Control",
                                           "Q" = "Least", "R" = "Most",
                                           "S" = "Most", "T" = "Control", 
                                           "U" = "Least", "V" = "Most", 
                                           "W" = "Most", "X" = "Most",
                                           "Y" = "Least", "Z" = "Least"))

## Calculate the time when each fish was infected

# Add the parasite intensity data to the metric dataframe
node_metrics$time <- as.integer(node_metrics$time)
inf_node_metrics <- left_join(node_metrics, infections, by = c("shoal",
                                                               "time",
                                                               "fish_id"))

# Identify the fish that are infected at each time step
inf_node_metrics$fish_infection_status <- cut(inf_node_metrics$parasite_intensity, c(-Inf,0,Inf))
levels(inf_node_metrics$fish_infection_status) <- c("Non-infected", "Infected")

# Identify the fish that will be infected at some point in the time series
infected_ever <- unique(subset(inf_node_metrics, fish_infection_status == "Infected")[,c(5:7)])

# Create a factor that indicates whether a fish ever gets infected in the experiment
inf_node_metrics$fish_infection_status_ever <- rep("Non-infected", nrow(inf_node_metrics))

for (i in 1:nrow(infected_ever)){
  fish <- infected_ever[i,]
  index <- rownames(inf_node_metrics[inf_node_metrics$shoal == fish$shoal & inf_node_metrics$fish_id == fish$fish_id,])
  inf_node_metrics[index,"fish_infection_status_ever"] <- "Infected"
}

# Make sure those fish where data doesn't exist on some time steps are NAs
inf_node_metrics[which(is.na(inf_node_metrics$fish_infection_status),
                       arr.ind=TRUE),"fish_infection_status_ever"] <- NA


#### Motif frequencies #### 

# Run the extract motifs functions (see function for description of methods)
motif_metrics <- data.frame()
for (i in 1:length(edges_list)) { 
  
  # Run the function and extract the counts (the 2nd object in the list)
  motif_counts <- motif.detector(edges_list[[i]])[[2]]
  
  # Extract the information on the shoal and time step for each network
  motif_metrics[i,"shoal"] <- str_split_fixed(names(edges_list), 
                                              pattern = "\\.", n = 2)[i,1]
  motif_metrics[i,"time"] <- str_split_fixed(names(edges_list), 
                                             pattern = "\\.", n = 2)[i,2]
  
  # Collate the data from the object
  motif_metrics[i,"os_n"] <- motif_counts[,"os_n"]
  motif_metrics[i,"is_n"] <- motif_counts[,"is_n"]
  motif_metrics[i,"inf_n"] <- motif_counts[,"inf_n"]
  motif_metrics[i,"inf_q"] <- motif_counts[,"inf_q"]
  motif_metrics[i,"os_q"] <- motif_counts[,"os_q"]
  motif_metrics[i,"is_q"] <- motif_counts[,"is_q"]
  
}

# Add before and after categories (before and after infections)
motif_metrics$infection <-  cut(as.numeric(motif_metrics$time), c(0,5,10))
levels(motif_metrics$infection) <- c("Before", "After")

# Add the treatments applied (infecting the most and least connected and a control)
motif_metrics$treatment <-  revalue(motif_metrics$shoal, 
                                    c("Aa" = "Least", "G" = "Least", 
                                      "H" = "Control", "I" = "Most",
                                      "J" = "Least", "K" = "Control", 
                                      "M" = "Most", "N" = "Least", 
                                      "O" = "Control", "P" = "Control",
                                      "Q" = "Least", "R" = "Most",
                                      "S" = "Most", "T" = "Control", 
                                      "U" = "Least", "V" = "Most", 
                                      "W" = "Most", "X" = "Most",
                                      "Y" = "Least", "Z" = "Least"))


#### Network metrics #### 

## Topological properties

# Calculate connectance
C <- unlist(lapply(nets_list, edge_density))

# Calculate reciprocity
reciprocity <- unlist(lapply(nets_list, reciprocity))

## Ratios of interactions between fish that are and are not infected

# Look at the ratio of interactions from high-low infected and low-high infected 
# Greater high-low interactions indicates shedding, greater low-high interactions indicates sharing
ratio_n <- NULL; ratio_q <- NULL; h2l_n <- NULL; l2h_n <- NULL; h2l_q <- NULL;
l2h_q <- NULL
for (i in names(edges_list)){ 
  
  edgelist <- edges_list[[i]] # Edgelist for each network
  
  sub <- edgelist[rowSums(edgelist[,c("lower_infection", "upper_infection")]) > 0,]

  if(nrow(sub) <= 1){ # If there is only one row (i.e., can't calculate a mean) or no data for infections
    ratio_n[[i]] <- NaN
    ratio_q[[i]] <- NaN
  } else { # Otherwise calculate the ratio of interactions for high-low and low-high interactions 
    h2l_n[[i]] <- nrow(sub[sub$lower_infection > sub$upper_infection,])
    l2h_n[[i]] <- nrow(sub[sub$lower_infection < sub$upper_infection,])
    h2l_q[[i]] <- sum(sub[sub$lower_infection > sub$upper_infection, "weight"])
    l2h_q[[i]] <- sum(sub[sub$lower_infection < sub$upper_infection, "weight"])
    
    # A high ratio would indicate shedding (i.e., more infected fish are trying to ditch parasites)
    # A low ratio would indicate sharing (i.e., less infect fish are trying to get parasites)
    ratio_n[[i]] <- h2l_n[[i]]/l2h_n[[i]] # With binary data for interactions
    ratio_q[[i]] <- h2l_q[[i]]/l2h_q[[i]] # With the edge weights for interactions
    }
}

# Bind together the different numbers of interactions (h2l and l2h)
directions <- data.frame(hightolow_n = unlist(h2l_n), 
                         lowtohigh_n = unlist(l2h_n), 
                         hightolow_q = unlist(h2l_q), 
                         lowtohigh_q = unlist(l2h_q))

# Extract information from names of rows 
directions$shoal <- unlist(map(str_split(rownames(directions), "\\."), 1))
directions$time <- as.integer(unlist(map(str_split(rownames(directions), "\\."), 2)))

# Create ratio objects
directionality_n <- unlist(ratio_n); directionality_q <- unlist(ratio_q)

# Inf values are produced when there are no lo2hi interactions or -Inf for the 
# opposite (so we will give them an arbitrary value above the maximum
directionality_n[is.infinite(directionality_n)] <- NA
directionality_q[is.infinite(directionality_n)] <- NA

# Make dataframe for network metrics
net_metrics <- data.frame(C, reciprocity, directionality_n, directionality_q)

# Extract information from names of rows 
net_metrics$shoal <- unlist(map(str_split(rownames(net_metrics), "\\."), 1))
net_metrics$time <- unlist(map(str_split(rownames(net_metrics), "\\."), 2))

# Add before and after categories (before and after infections)
net_metrics$infection <-  cut(as.numeric(net_metrics$time), c(0,5,10))
levels(net_metrics$infection) <- c("Before", "After")

# Add the treatments applied (infecting the most and least connected and a control)
net_metrics$treatment <-  revalue(net_metrics$shoal, 
                                         c("Aa" = "Least", "G" = "Least", 
                                           "H" = "Control", "I" = "Most",
                                           "J" = "Least", "K" = "Control", 
                                           "M" = "Most", "N" = "Least", 
                                           "O" = "Control", "P" = "Control",
                                           "Q" = "Least", "R" = "Most",
                                           "S" = "Most", "T" = "Control", 
                                           "U" = "Least", "V" = "Most", 
                                           "W" = "Most", "X" = "Most",
                                           "Y" = "Least", "Z" = "Least"))


#### Getting data frames with complete data for analysis and plotting 

# Count the number of infected fish
infection_no <- aggregate(parasite_intensity ~ shoal + time, 
                          FUN = function(x){length(unique(x[x != 0]))}, 
                          data = infections)
colnames(infection_no) <- c("shoal", "time", "infected_fish_no")

# Convert motif data to long format for plotting
motifs_long <- melt(dplyr::select(motif_metrics, -is_n, -is_q),
                    id.vars = c("shoal", "time", "infection", "treatment"))
motifs_plotting <- subset(motifs_long, infection == "After")
motifs_plotting$time <- as.numeric(motifs_plotting$time)
motifs_plotting$variable <- ordered(motifs_plotting$variable,
                                    levels = c("inf_n", "inf_q", "os_n",
                                               "os_q"))

# Add the number of infected species to the motif metrics dataframe 
motif_wide <- left_join(motifs_plotting, infection_no, by = c("shoal", "time"))

# Get mean values of parasite intensity for shoals
infections_population <- aggregate(parasite_intensity ~ shoal + time, 
                                   FUN = "mean", data = infections)

# Merge the two data frames
inf_motifs_long <- left_join(motif_wide, infections_population, 
                             by = c("shoal", "time"))

# Merge the two data frames
inf_directions <- left_join(directions, infections_population, 
                             by = c("shoal", "time"))

# Merge the two data frames 
net_metrics$time <- as.integer(net_metrics$time)
inf_net_metrics <- left_join(net_metrics, infections_population, 
                             by = c("shoal", "time"))


## Remove unneccessary objects 

# Remove objects
rm(betweenness, closeness, C, directionality_n, directionality_q, i, index,
   reciprocity, windegree, woutdegree, motif.detector, targets, sub, ratio_n,
   ratio_q, fish, motif_counts, infected_ever, edges, edgelist, 
   edges_list, l2h_n, h2l_n, l2h_q, h2l_q, nets_I_a, nets_I_b, nets_I_ba)
