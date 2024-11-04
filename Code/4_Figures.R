# Social network dynamics in response to parasitism (guppy-gyro model) #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line # 


## 3 - Statistical analysis



#### Setup ####

## Get the environment ready for loading the data

# Clear environment
rm(list=ls())

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Social networks in guppies (Scientific Reports)")

# Load the previous scripts 
source("Code/4_Statistical-analysis.R")


#### Figure 1 #### 

## Plot the out degree for the different fish in the experiments

# ggplot rendering for weighted in degree
plot1a <- ggplot(aes(x = infection, y = windegree, fill = treatment),
                data = node_metrics) + 
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.5, position = position_dodge(0.5), inherit.aes = F, 
               data = node_metrics) +
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = node_metrics) + 
  theme_bw() + 
  theme(legend.position = c(0.36, 0.75), legend.margin=margin(c(1,1,1,1)),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected")) + 
  ylab("Weighted in degree") + 
  xlab("") + 
  ggtitle("a")
plot1a

# ggplot rendering for weighted out degree
plot1b <- ggplot(aes(x = infection, y = woutdegree, fill = treatment),
                data = node_metrics) + 
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.5, position = position_dodge(0.5), inherit.aes = F, 
               data = node_metrics) +
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = node_metrics) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected")) + 
  ylab("Weighted out degree") + 
  xlab("") + 
  ggtitle("b")
plot1b

# ggplot rendering for betweenness
plot1c <- ggplot(aes(x = infection, y = betweenness, fill = treatment),
                data = node_metrics) + 
  geom_violin(position = position_dodge(0.5)) +
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.5, position = position_dodge(0.5), inherit.aes = F, 
               data = node_metrics) +
  stat_summary(aes(x = infection, y = woutdegree, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = node_metrics) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "black"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected"), 
                       option = "plasma") + 
  ylab("Betweenness") + 
  xlab("") + 
  ggtitle("c")
plot1c

# ggplot rendering for closeness
plot1d <- ggplot(aes(x = infection, y = closeness, fill = treatment),
                 data = node_metrics) + 
  geom_violin(position = position_dodge(0.5)) +
  stat_summary(aes(x = infection, y = closeness, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.5, position = position_dodge(0.5), inherit.aes = F, 
               data = node_metrics) +
  stat_summary(aes(x = infection, y = closeness, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = node_metrics) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "black"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected")) + 
  ylab("Closeness") + 
  xlab("") + 
  ggtitle("d")
plot1d

# Plot the change in weighted out degree to compare between fish that eventually become infected and those that don't 
inf_node_metrics$fish_infection_status_ever <- as.factor(inf_node_metrics$fish_infection_status_ever)
inf_node_metrics$fish_infection_status_ever <- relevel(inf_node_metrics$fish_infection_status_ever, "Non-infected")
plot1e <- ggplot(aes(x = fish_infection_status_ever, # Whether a fish goes on to become infected 
                     y = woutdegree, # Weight out degree
                     fill = infection), # Before or after the infection (i.e., day 1-5 or 6-10)
                # Remove some fish where parasite load wasn't counted
               data = subset(inf_node_metrics, 
                             !is.na(fish_infection_status_ever))) + 
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = fish_infection_status_ever, y = woutdegree, 
                   group = infection), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.2, position = position_dodge(0.5), inherit.aes = F, 
               data = inf_node_metrics) +
  stat_summary(aes(x = fish_infection_status_ever, y = woutdegree, 
                   group = infection), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = inf_node_metrics) + 
  scale_fill_manual(values = c("grey50", "darkred"), 
                    labels = c("Before", "After"), name = "Time period") + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  ylab("Weighted out degree") + 
  xlab("") + 
  ggtitle("f")
plot1e

# Plot the change in weighted in degree to compare between fish that eventually
# become infected and those that don't 
plot1f <- ggplot(aes(x = fish_infection_status_ever, y = windegree, 
                     fill = infection),
                 data = subset(inf_node_metrics, 
                               !is.na(fish_infection_status_ever))) +
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = fish_infection_status_ever, y = windegree, 
                   group = infection), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.2, position = position_dodge(0.5), inherit.aes = F, 
               data = inf_node_metrics) +
  stat_summary(aes(x = fish_infection_status_ever, y = windegree, 
                   group = infection), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = inf_node_metrics) + 
  scale_fill_manual(values = c("grey50", "darkred"), 
                    labels = c("Before", "After"), name = "Time period") + 
  theme_bw() + 
  theme(legend.position = c(0.755,0.8), legend.margin=margin(c(1,1,1,1)),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  ylab("Weighted in degree") + 
  xlab("") + 
  ggtitle("e")
plot1f

# Plot parasite intensity versus time
plot1g <- ggplot(aes(y = parasite_intensity, x = time), 
                 data = subset(inf_node_metrics, 
                        !is.na(fish_infection_status_ever))) +  
  geom_point(aes(fill = fish_infection_status), size = 2, pch = 21) + 
  scale_fill_manual(values = c("white", "black"), 
                    name = "Fish infection status") + 
  theme_bw() + 
  theme(legend.position = c(0.17, 0.68), legend.margin=margin(c(1,1,1,1)),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(1,10,1)) + 
  ylab("Parasite intensity (n)") + 
  xlab("Time (days)") + 
  ggtitle("g")
plot1g

# Plot in versus out degree for non infected fish 
plot1h <- ggplot(aes(x = windegree, y = woutdegree, 
                     fill = fish_infection_status),
                 data = subset(inf_node_metrics, 
                               fish_infection_status == "Non-infected")) +
  geom_point(fill = "white", pch = 21, size = 2) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        legend.direction = "horizontal") +
  ylab("Weighted out degree") + 
  xlab("Weighted in degree") + 
  ylim(0,25) + 
  xlim(0,18) + 
  ggtitle("h")
plot1h

# Plot in versus out degree for infected fish 
plot1i <- ggplot(aes(x = windegree, y = woutdegree, 
                     fill = fish_infection_status),
                 size = 2, pch = 21,
                 data = subset(inf_node_metrics, 
                               fish_infection_status == "Infected")) +
  geom_point(fill = "black", pch = 21, size = 2) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        legend.direction = "horizontal") +
  ylab("Weighted out degree") + 
  xlab("Weighted in degree") +
  ylim(0,25) + 
  xlim(0,18) + 
  ggtitle("i")
plot1i

# Arrange the plots (647 x 950)
grid.arrange(plot1a, plot1b, plot1c, plot1e, plot1f, plot1g, plot1h,
             plot1i, layout_matrix = rbind(c(1,2,3,NA),
                                           c(1,2,3,NA),
                                           c(1,2,3,NA),
                                           c(1,2,3,NA),
                                           c(5,4,NA,NA),
                                           c(5,4,NA,NA),
                                           c(5,4,NA,NA),
                                           c(5,4,NA,NA),
                                           c(6,6,NA,NA),
                                           c(6,6,NA,NA),
                                           c(6,6,NA,NA),
                                           c(7,8,NA,NA),
                                           c(7,8,NA,NA),
                                           c(7,8,NA,NA)))


#### Figure 2 ####

# Create labels 
motif_labels <- c(`inf_n` = "Asymmetric (n)", `inf_q` = "Asymmetric (q)",
                  `os_n` = "Out-star (n)", `os_q` = "Out-star (q)")

# Plot the relationship between parasite intensity and the motif counts
plot2a <- ggplot(aes(x = infected_fish_no, y = value, 
                     fill = as.factor(time)), data = inf_motifs_long) + 
  geom_point(pch = 21, size = 2) + 
  geom_smooth(aes(x = infected_fish_no, y = value), inherit.aes = F,
              colour = "black", method = "lm") + 
  theme_bw() +
  theme(legend.position = "NA",  
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 8)) + 
  xlab("Number of infected fish") + 
  ylab("Motif count") +
  scale_fill_viridis_d(name = "Time (day)", option = "magma") +
  facet_wrap(~variable, scales = "free_y", 
             labeller = labeller(variable = motif_labels)) + 
  ggtitle("a")
plot2a

# Plot the counts of the data
plot2b <- ggplot(aes(x = variable, y = value/infected_fish_no, 
                     fill = as.factor(time)), 
                 data = subset(motif_wide, treatment != "Control")) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  scale_fill_viridis_d(name = "Time (day)", option = "magma") +
  theme_bw() + 
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 8)) + 
  ylab("Motif count/ number of infected fish") + 
  xlab("Motif type") + 
  scale_x_discrete(labels=c("Asymmetric (n)", "Asymmetric (q)", "Out-star (n)", 
                            "Out-star (q)")) + 
  facet_wrap(~treatment, scales = "free_y") +
  ggtitle("b")
plot2b

grid.arrange(plot2a, plot2b, nrow = 2)


#### Figure 3 ####

## Plot the change in edge density and reciprocity 

# ggplot rendering
plot3a <- ggplot(aes(x = infection, y = C, fill = treatment),
                data = net_metrics) + 
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = infection, y = C, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.25, position = position_dodge(0.5), inherit.aes = F, 
               data = net_metrics) +
  stat_summary(aes(x = infection, y = C, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = net_metrics) + 
  theme_bw() + 
  theme(legend.position = c(0.36,0.78), 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected")) + 
  ylab("Edge density") + 
  xlab("Infection") + 
  ggtitle("a")
plot3a

# ggplot rendering
plot3b <- ggplot(aes(x = infection, y = reciprocity, fill = treatment),
                data = net_metrics) + 
  geom_violin(position = position_dodge(0.5)) + 
  stat_summary(aes(x = infection, y = C, group = treatment), 
               geom = "errorbar", fun.data = "mean_cl_boot",
               width = 0.25, position = position_dodge(0.5), inherit.aes = F, 
               data = net_metrics) +
  stat_summary(aes(x = infection, y = C, group = treatment), 
               geom = "point", fun = "mean", position = position_dodge(0.5), 
               size = 2, pch = 21, inherit.aes = F, fill = "white",
               data = net_metrics) + 
  theme_bw() + 
  theme(legend.position = "NA", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(name = "Treatment", labels = c("Control", 
                                                      "Least connected", 
                                                      "Most connected")) + 
  ylab("Reciprocity") + 
  xlab("Infection") + 
  ggtitle("b")
plot3b

# Before network for shoal I
net1 <- ggplot(data = fortify(nets_I_b)) + 
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, 
                 colour = as.numeric(weight)), 
              arrow = arrow(length = unit(1, "lines"),
                            type = "closed", angle = 25), curvature = 0.2) + 
  geom_nodes(aes(x = x, y = y, fill = as.character(infection)), size = 4, 
             pch = 21, show.legend = F) +
  theme_void() + 
  scale_fill_manual(labels = c("Uninfected"), values = c("darkgrey")) + 
  scale_colour_viridis_c(option = "plasma", name = "Weight") + 
  scale_linewidth_continuous(limits = c(1,20)) + 
  theme(legend.position = c(0.75,0.95), legend.direction = "horizontal", 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm")) + 
  ggtitle("c")
net1

# After network for shoal I
I_a <- left_join(fortify(nets_I_a), para_I_a, by = c("name" = "fish_id"))
I_a$para_int <- sqrt(I_a$parasite_intensity) + 4

net2 <- ggplot(data = I_a) + 
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend,
                 colour = as.numeric(weight)),
             arrow = arrow(length = unit(1, "lines"),
                           type = "closed", angle = 25), curvature = 0.2) + 
  geom_nodes(aes(x = x, y = y, fill = as.character(infection.x), 
                 size = para_int), pch = 21) +
  geom_node_label(aes(label = parasite_intensity), nudge_y = 0.08) + 
  theme_void() + 
  scale_fill_manual(labels = c("Uninfected", "Infected"), 
                    values = c("darkgrey", "darkred")) + 
  scale_colour_viridis_c(option = "plasma", name = "Weight") + 
  scale_size_continuous(limits = c(1,16)) + 
  theme(legend.position = "NA", 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm")) + 
  ggtitle("d")
net2

# ggplot rendering
plot3e <- ggplot(aes(x = lowtohigh_n, y = hightolow_n, fill = as.factor(time)),
                 data = directions) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point(pch = 21, size = 3) + 
  theme_bw() + 
  theme(legend.position = "NA",
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(option = "magma", name = "Day") + 
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) + 
  ylab("High-low infection interactions (n)") + 
  xlab("Low-high infection interactions (n)") + 
  ggtitle("e")
plot3e

# ggplot rendering
plot3f <- ggplot(aes(x = lowtohigh_q, y = hightolow_q, fill = as.factor(time)),
                 data = directions) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  geom_point(pch = 21, size = 3) + 
  theme_bw() + 
  theme(legend.position = "right", 
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        legend.background = element_rect(colour = "white"), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) + 
  scale_fill_viridis_d(option = "magma", name = "Day") + 
  ylab("High-low infection interactions (q)") + 
  xlab("Low-high infection interactions (q)") + 
  ggtitle("f")
plot3f

plot3e <- ggplotGrob(plot3e); plot3f <- ggplotGrob(plot3f)
plot3e$widths <- plot3f$widths

# Multiplot
grid.arrange(plot3a, plot3b, net1, net2, plot3e, plot3f, 
             layout_matrix = rbind(c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA,NA), 
                                   c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA,NA),
                                   c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA,NA),
                                   c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA,NA),
                                   c(3,3,3,3,3,3,3,4,4,4,4,4,4,4),
                                   c(3,3,3,3,3,3,3,4,4,4,4,4,4,4),
                                   c(3,3,3,3,3,3,3,4,4,4,4,4,4,4),
                                   c(5,5,5,5,5,5,5,6,6,6,6,6,6,6),
                                   c(5,5,5,5,5,5,5,6,6,6,6,6,6,6),
                                   c(5,5,5,5,5,5,5,6,6,6,6,6,6,6))) 

