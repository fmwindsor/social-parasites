# Social network dynamics in response to parasitism (guppy-gyro model) #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## Motif detector function


# A function to detect motifs of interest including infected fish, adapted from:

# Tavella et al. (2022) Using motifs in ecological networks to identify the role
# of plants in crop margins for multiple agriculture functions. Agriculture,
# Ecosystems & Environment, 331, 107912.


# Function definition
motif.detector <- function(edgelist){
	# Input = edgelist containing interactions
	# Output = a list of motifs

  # Subset the edgelists to only those direct interactions with infected fish
  infected_edges <- edgelist[edgelist$lower_infection > 0 |
                               edgelist$upper_infection > 0,] 
  
  # Get lists of infected fish (lower nodes)
  infected_lower <- unique(infected_edges[infected_edges$lower_infection > 0,
                                          "lower"])
  
  # Get lists of infected fish (higher nodes)
  infected_upper <- unique(infected_edges[infected_edges$upper_infection > 0,
                                          "upper"])
  
  if(nrow(infected_edges) >= 1){
    
    
    ## Dyad motif (interaction between one infected and one uninfected)

    # Store the edges for the individual interactions
    inf_dframe <- infected_edges
  
    # Create null objects to store the count and weighted count data
    n_inf_motif <- NULL; q_inf_motif <- NULL
    
    # Calculate the number of motifs as well as a sum of their weights
    for (i in unique(inf_dframe$lower)){ 
      n_inf_motif[i] <- nrow(inf_dframe[inf_dframe$lower == i,])
      q_inf_motif[i] <- sum(inf_dframe[inf_dframe$lower == i, "weight"])
    }
    
    # Calculate the count and weighted counts 
    inf_count <- sum(n_inf_motif); inf_value <- sum(q_inf_motif)
 
    
    ## Triad motif 1 (interactions between one infected and two other fish)
    
    if(length(infected_lower) >= 1){ 

      # Create a null object to store the data
      os_motif <- NULL 
    
      for (i in infected_lower){
        # If there are more than two interactions starting from one infected fish
        if(nrow(infected_edges[infected_edges$lower == i,]) >= 2){ 
          # Multiple fish interacting with the infected fish
          os_motif[[i]] <- infected_edges[infected_edges$lower == i,] 
          }
      }
      # If there is not a motif detected
      if(is.null(os_motif)){
        # Create an empty/null object
        os_dframe <- NULL
        } else{
          # Bind the lists of interactions together
          os_dframe <- as.data.frame(do.call(rbind, os_motif))
        }
  
      # Count the occurrences of the motifs (n = unweighted, q = weighted)
      n_os_motif <- NULL; q_os_motif <- NULL; 
      
      for (j in unique(os_dframe$lower)){ 
        n_os_motif[j] <-  nrow(os_dframe[os_dframe$lower == j,])
        q_os_motif[j] <- sum(os_dframe[os_dframe$lower == j, "weight"])
      }
      os_count <- sum(n_os_motif)
      os_value <- sum(q_os_motif)
    
      } else{ 
        os_dframe <- NULL; os_count <- 0; os_value <- 0
      }
    
    
    ## Triad motif 2 (interactions between two infected fish and one other fish)
    
    if(length(infected_lower) >= 2){ 
   
      # Create a null object to store the data 
      is_motif <- NULL
      
      # For every upper node (i.e., the potential for an instar)
      for (k in infected_upper){
        # If there are more than two species interacting with the upper node
        if(nrow(infected_edges[infected_edges$upper == k,]) >= 2){
          # Fish interacting with the infected fish
          is_motif[[k]] <- infected_edges[infected_edges$upper == k,] 
        }
      }
      
      if(is.null(is_motif)){ 
        is_dframe <- NULL
        } else{
          # Bind the lists of interactions together
          is_dframe <- as.data.frame(do.call(rbind, is_motif))
        }
      
      # Create null objects to store data
      n_is_motif <- NULL; q_is_motif <- NULL
      
      # Count the occurrences of the in-star motifs
      for (i in unique(is_motif$upper)){
        n_is_motif[i] <- nrow(is_dframe[is_dframe$upper == i])
        q_is_motif[i] <- sum(is_dframe[is_dframe$upper == i, "weight"])
      }
      is_count <- sum(n_is_motif); is_value <- sum(q_is_motif)
      
      } else{
        is_dframe <- NULL; is_count <- 0; is_value <- 0
      }
    
    # Collate the count data into a dataframe
    motif_summary <- data.frame(os_n = os_count, os_q = os_value,
                                inf_n = inf_count, inf_q = inf_value,
                                is_n = is_count, is_q = is_value)
    
    # Collate and return the motif edgelists and count data
    return(list(infected_edges, motif_summary))
    
    } else{ # If there are no motifs print null
      return(list(NULL, data.frame(os_n = 0, is_n = 0, inf_n = 0, inf_q = 0,
                                   os_q = 0, is_q = 0)))
    }
}

