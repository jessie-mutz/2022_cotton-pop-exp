## USDA cotton experiment
## Figures
## Code written and annotated by JM (last updated 2024-11-08)

require(ggplot2)
require(gridExtra)
require(ggtext)

## =============================================================================
## Figure 1. Plant responses with populations means (across timepoints and herbivore densities)
## =============================================================================

# Read in data for defenses and resistance.
full_data <- read.csv("cleaned-data/full-data.csv")

# Add a column for population ID (Cage_pos) as integers 1-70.
full_data$Pop_ID <- with(full_data, substr(Plant_ID - 1, start = 2, stop = 3))
full_data$Pop_ID <- with(full_data, as.numeric(Pop_ID) + 1)

# Scale leaf size measurements (initial leaf area or midvein length) for use as
#  covariates
full_data$Initial_leaf_area_sc <- scale(full_data$Initial_leaf_area)[,1]
full_data$Midvein_length_sc <- scale(full_data$Midvein_length)[,1]

# Read in summaries from single trait defense or resistance models.
rgr_summary <- read.csv("output/full-rgr-summary-min5.csv")
con_summary <- read.csv("output/full-con-summary-min5.csv")
tri_summary <- read.csv("output/full-tri-summary-min5.csv")
H10_summary <- read.csv("output/full-H10-summary-min5.csv")
H36_summary <- read.csv("output/full-H36-summary-min5.csv")
gna_summary <- read.csv("output/full-gland-area-summary-min5.csv")

## Standardize observations to mean leaf size

# For each response, determine whether there is an effect of leaf size
rgr_summary[which(rgr_summary$X == "slope_resp"), ]
con_summary[which(con_summary$X == "slope_resp"), ] #no effect here
tri_summary[which(tri_summary$X == "slope_resp"), ]
H10_summary[which(H10_summary$X == "slope_resp"), ] #no effect here
H36_summary[which(H36_summary$X == "slope_resp"), ]
#gnd_summary[which(gnd_summary$X == "slope_resp"), ]
gna_summary[which(gna_summary$X == "slope_resp"), ]

# For response variables with significant effect of leaf size, standardize trait 
#  values to mean leaf size
full_data$rgr_sizestd <- with(full_data, RGR_mean - 
                                rgr_summary[which(rgr_summary$X == "slope_resp"), 'mean']*Initial_leaf_area_sc)
full_data$con_sizestd <- with(full_data, ifelse(Mean_area_consumed == 0, NA,
                                                log(Mean_area_consumed)))
full_data$Trichomes_sizestd <- with(full_data, sqrt(Trichomes_std_cm2) - 
                                      tri_summary[which(tri_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$H_10.5_sizestd <- with(full_data, log(H_10.5_relabun))
full_data$H_3.6_sizestd <- with(full_data, log(H_3.6_relabun) - 
                                  H36_summary[which(H36_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$Gland_prop_sizestd <- with(full_data, sqrt(Gland_prop) - 
                                       gna_summary[which(gna_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)

# Read in population means for each response.
rgr_pop <- read.csv("output/full-rgr-pop-means-min5.csv")
con_pop <- read.csv("output/full-con-pop-means-min5.csv")
tri_pop <- read.csv("output/full-tri-pop-means-min5.csv")
gnd_pop <- read.csv("output/full-gland-dens-pop-means-min5.csv")
gna_pop <- read.csv("output/full-gland-area-pop-means-min5.csv")
H10_pop <- read.csv("output/full-H10-pop-means-min5.csv")
H36_pop <- read.csv("output/full-H36-pop-means-min5.csv")

## Function for plotting size standardized plant level responses with population
##  means (+/- 89% CI). Response ~ days (x axis) x herbivore density (dodged points).
plot1 <- function(pop_df, model_df, response) {
  ggplot(data = pop_df, aes(x = Day, y = mean, 
                            group = as.factor(Herb_add),
                            fill = as.factor(Herb_add),
                            shape = as.factor(Herb_add))) +
    geom_point(data = model_df, aes(x = model_df[,'Day'], y = model_df[,response], 
                                    color = as.factor(Herb_add)),
               position = position_dodge(width = 1.5), size = 1.5, alpha = 0.35,
               stroke = 0, color = "transparent") +
    geom_errorbar(aes(x = Day, ymin = CI_low_95, ymax = CI_high_95, 
                      color = as.factor(Herb_add)),
                  position = position_dodge(width = 1.5), width = 1, linewidth = 0.3) +
    geom_point(position = position_dodge(width = 1.5), size = 2.1, 
               color = "black", stroke = 0.1) +
    theme_bw() + theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      axis.title.x = ggtext::element_markdown(size = 12.5),
      axis.title.y = ggtext::element_markdown(size = 12.5),
      axis.text = element_text(size = 12, color = "black"),
      legend.key.height = unit(0.7, 'cm'),
      legend.key.width = unit(0.7, 'cm'),
      legend.text = element_text(size = 12),
      legend.title = ggtext::element_markdown(size = 12, hjust = 0.5),
      legend.margin = margin(2,5,2,2),
      plot.tag.location = "panel",
      plot.tag = element_text(margin = margin(0,0,0,4)),
      plot.margin = margin(10,10,10,12),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
    ) +
    labs(x = "Days following introduction of<br>experimental *S. exigua*",
         y = "sqrt Trichomes cm<sup>-2</sup>", 
         fill = "Experimental<br>*S. exigua*<br>density<br>(per plant)", 
         shape = "Experimental<br>*S. exigua*<br>density<br>(per plant)",
         color = "Experimental<br>*S. exigua*<br>density<br>(per plant)") +
    guides(color = guide_legend(override.aes = list(size = 3.5))) +
    scale_x_continuous(breaks = c(0,2,4,7,12,21), limits = c(0, 22), expand = c(0,0)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_shape_manual(values = c(23, 22, 24, 21))
}

## Trichome density as response.
# Set up dataframe.
tri_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'Trichomes_sizestd', 'Midvein_length_sc')]
tri_sub <- na.omit(tri_sub)
# Include only populations with at least 5 plant observations.
tri_sub <- tri_sub[which(tri_sub$Pop_ID %in% names(which(table(tri_sub$Pop_ID) >= 5))),]

tri_pred <- read.csv("output/full-tri-summary-min5.csv")
tri_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(tri_pred[1:5, 'mean'], tri_pred[1:5, 'mean'] + 
             tri_pred[grepl('pred', tri_pred$X), 'mean'])
)

# Using plotting function.
plot1(tri_pop, tri_sub, 'Trichomes_sizestd') +
  geom_line(data = tri_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")

# Use plotting function.
plot1(tri_pop, tri_sub, 'Trichomes_sizestd') 
  
# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/trichomes.pdf", width = 0.9*6.09, height = 0.9*4.21)
plot1(tri_pop, tri_sub, 'Trichomes_sizestd')
dev.off()
png(file = "figures/trichomes.png", width = 0.9*6.09, height = 0.9*4.21, 
    units = "in", res = 600)
plot1(tri_pop, tri_df, 'Trichomes_sizestd')
dev.off()


## Hemigossypolone 10.5 relative abundance as response.
# Set up dataframe.
H10_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'H_10.5_sizestd', 'Midvein_length_sc')]
H10_sub <- na.omit(H10_sub)
# Include only populations with at least 5 plant observations.
H10_sub <- H10_sub[which(H10_sub$Pop_ID %in% names(which(table(H10_sub$Pop_ID) >= 5))),]

H10_pred <- read.csv("output/full-H10-summary-min5.csv")
H10_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(H10_pred[1:5, 'mean'], H10_pred[1:5, 'mean'] + 
             H10_pred[grepl('pred', H10_pred$X), 'mean'])
)

# Using plotting function.
plot1(H10_pop, H10_sub, 'H_10.5_sizestd') +
  labs(y = "ln HGQ 10.5")  +
  geom_line(data = H10_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/H10.5-plot.pdf", width = 0.9*6.09, height = 0.9*4.21)
#plot1(H10_pop, H10_sub, 'H_10.5_sizestd') +
#  labs(y = "ln HGQ 10.5")
dev.off()
png(file = "figures/H10.5-plot.png", width = 0.9*6.09, height = 0.9*4.21,
    units = "in", res = 600)
#plot1(H10_pop, H10_sub, 'H_10.5_sizestd') +
#  labs(y = "ln HGQ 10.5")
dev.off()



## Hemigossypolone 3.6 relative abundance as response.
# Set up dataframe.
H36_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'H_3.6_sizestd', 'Midvein_length_sc')]
H36_sub <- na.omit(H36_sub)
# Include only populations with at least 5 plant observations.
H36_sub <- H36_sub[which(H36_sub$Pop_ID %in% names(which(table(H36_sub$Pop_ID) >= 5))),]

H36_pred <- read.csv("output/full-H36-summary-min5.csv")
H36_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(H36_pred[1:5, 'mean'], H36_pred[1:5, 'mean'] + 
             H36_pred[grepl('pred', H36_pred$X), 'mean'])
)

# Using plotting function.
plot1(H36_pop, H36_sub, 'H_3.6_sizestd') +
  labs(y = "ln HGQ 3.6")  +
  geom_line(data = H36_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/H3.6-plot.pdf", width = 0.9*6.09, height = 0.9*4.21)
#plot1(H10_pop, H36_sub, 'H_3.6_sizestd') +
#  labs(y = "ln HGQ 3.6")
dev.off()
png(file = "figures/H10.5-plot.png", width = 0.9*6.09, height = 0.9*4.21,
    units = "in", res = 600)
#plot1(H10_pop, H36_sub, 'H_3.6_sizestd') +
#  labs(y = "ln HGQ 3.6")
dev.off()


## Gland proportional leaf area as response.
# Set up dataframe.
gna_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'Gland_prop_sizestd', 'Midvein_length_sc')]
gna_sub <- na.omit(gna_sub)
# Include only populations with at least 5 plant observations.
gna_sub <- gna_sub[which(gna_sub$Pop_ID %in% names(which(table(gna_sub$Pop_ID) >= 5))),]


gna_pred <- read.csv("output/full-gland-area-summary-min5.csv")
gna_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(gna_pred[1:5, 'mean'], gna_pred[1:5, 'mean'] + 
             gna_pred[grepl('pred', gna_pred$X), 'mean'])
)

# Use plotting function.
plot1(gna_pop, gna_sub, 'Gland_prop_sizestd') +
  labs(y = "sqrt Prop. gland area") +
  geom_line(data = gna_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/gland-area-plot.pdf", width = 0.9*6.09, height = 0.9*4.21)
#plot1(gna_pop, gna_sub, 'Gland_prop_sizestd') +
#  labs(y = "sqrt Prop. gland area")
dev.off()
png(file = "figures/gland-area-plot.png", width = 0.9*6.09, height = 0.9*4.21,
    units = "in", res = 600)
#plot1(gna_pop, gna_sub, 'Gland_prop_sizestd') +
#  labs(y = "sqrt Prop. gland area")
dev.off()


## Leaf area consumed as response.
# Set up dataframe.
con_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'con_sizestd', 'Initial_leaf_area_sc')]
con_sub <- na.omit(con_sub)
# Include only populations with at least 5 plant observations.
con_sub <- con_sub[which(con_sub$Pop_ID %in% names(which(table(con_sub$Pop_ID) >= 5))),]

con_pred <- read.csv("output/full-con-summary-min5.csv")
con_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(con_pred[1:5, 'mean'], con_pred[1:5, 'mean'] + 
             con_pred[grepl('pred', con_pred$X), 'mean'])
)

# Use plotting function.
plot1(con_pop, con_sub, 'con_sizestd') +
  labs(y = "ln Bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)") +
  geom_line(data = con_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/area-consumed-plot.pdf", width = 0.9*7.09, height = 0.9*4.21)

dev.off()
png(file = "figures/area-consumed-plot.png", width = 0.9*7.09, height = 0.9*4.21,
    units = "in", res = 600)

dev.off()

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/area-consumed-plot-no-lines.pdf", width = 0.9*7.09, height = 0.9*4.21)
plot1(con_pop, con_sub, 'con_sizestd') +
  labs(y = "ln Leaf area consumed (cm<sup>2</sup>)") 
dev.off()
png(file = "figures/area-consumed-plot-no-lines.png", width = 0.9*7.09, height = 0.9*4.21,
    units = "in", res = 600)
plot1(con_pop, con_sub, 'con_sizestd') +
  labs(y = "ln Leaf area consumed (cm<sup>2</sup>)") 
dev.off()


## Bioassay caterpillar RGR as response.
# Set up dataframe.
rgr_sub <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day',
                       'rgr_sizestd', 'Initial_leaf_area_sc')]
rgr_sub <- na.omit(rgr_sub)
# Include only populations with at least 5 plant observations.
rgr_sub <- rgr_sub[which(rgr_sub$Pop_ID %in% names(which(table(rgr_sub$Pop_ID) >= 5))),]

rgr_pred <- read.csv("output/full-rgr-summary-min5.csv")
rgr_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(rgr_pred[1:5, 'mean'], rgr_pred[1:5, 'mean'] + 
             rgr_pred[grepl('pred', rgr_pred$X), 'mean'])
)

# Using plotting function.
plot1(rgr_pop, rgr_sub, 'rgr_sizestd') +
  labs(y = "Bioassay caterpillar RGR (g g<sup>-1</sup> d<sup>-1</sup>)")

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/RGR-plot.pdf", width = 6.09, height = 4.21)
plot1(rgr_pop, rgr_sub, 'rgr_sizestd') +
  labs(y = "Caterpillar RGR (g g<sup>-1</sup> d<sup>-1</sup>)")
dev.off()
png(file = "figures/RGR-plot.png", width = 6.09, height = 4.21,
    units = "in", res = 600)
plot1(rgr_pop, rgr_sub, 'rgr_sizestd') +
  labs(y = "Caterpillar RGR (g g<sup>-1</sup> d<sup>-1</sup>)")
dev.off()


  
  
  tri_plot <-
    plot1(tri_pop, tri_sub, 'Trichomes_sizestd') +
    theme(legend.position = "none") +
    labs(y = "sqrt Trichomes cm<sup>-2</sup>", tag = "A") +
    annotate("text", x = c(21), y = c(Inf), vjust = 1,
             label = c("+")) +
  geom_line(data = tri_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  H10_plot <- 
    plot1(H10_pop, H10_sub, 'H_10.5_sizestd') +
    theme(legend.position = "none") +
    labs(y = "ln HGQ 10.5", tag = "D") +
    annotate("text", x = c(4,7,12), y = c(Inf, Inf, Inf), vjust = 1,
             label = c("*", "*", "*")) +
  geom_line(data = H10_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  H36_plot <- 
    plot1(H36_pop, H36_sub, 'H_3.6_sizestd') +
    theme(legend.position = "none") +
    labs(y = "ln HGQ 3.6", tag = "C")  +
    annotate("text", x = c(4,12), y = c(Inf, Inf), vjust = 1, label = c("*", "+")) +
  geom_line(data = H36_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  
  gna_plot <-
    plot1(gna_pop, gna_sub, 'Gland_prop_sizestd') +
    theme(legend.position = "none") +
    labs(y = "sqrt Prop. gland area", tag = "B") +
    annotate("text", x = c(7,21), y = c(Inf, Inf), vjust = 1, label = c("*", "*")) +
  geom_line(data = gna_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  con_plot <-
    plot1(con_pop, con_sub, 'con_sizestd') +
    theme(legend.position = "none")+
    labs(y = "ln Bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)", tag = "E") +
    annotate("text", x = c(2,4,7,12), y = c(Inf, Inf, Inf, Inf), vjust = 1,
             label = c("*", "*", "*", "*")) +
  geom_line(data = con_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  rgr_plot <-
    plot1(rgr_pop, rgr_sub, 'rgr_sizestd') +
    theme(legend.position = "none")+
    labs(y = "Bioassay caterpillar<br>RGR (g g<sup>-1</sup> d<sup>-1</sup>)", tag = "F") +
    annotate("text", x = c(4,7,21), y = c(Inf, Inf, Inf), vjust = 1,
             label = c("*", "*", "+")) +
  geom_line(data = rgr_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), guide = "none")
  
  
ggpubr::as_ggplot(ggpubr::get_legend(plot1(tri_pop, tri_sub, 'Trichomes_sizestd')))  

grid::grid.newpage()
ga <- arrangeGrob(
  arrangeGrob(
    grid::grid.text(label = "Defenses", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    grid::grid.text(label = "", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    tri_plot, gna_plot, H36_plot, H10_plot, nrow = 3, heights = c(0.1,1,1)),
  ggpubr::as_ggplot(ggpubr::get_legend(plot1(tri_pop, tri_sub, 'Trichomes_sizestd'))),
                    ncol = 2, widths = c(1,0.16) )
gb <- grid::grid.rect(width = 0.985, height = 0.985, 
                gp = grid::gpar(lwd = 1.5, col = "black", fill = NA))
#gc <- grid::grid.text(label = "Defenses", hjust = 0, vjust = 1, 
#                x = unit(0.02, "npc"), y = unit(0.98, "npc"), 
#                gp = grid::gpar(fontsize = 15, fontface = "bold"))
gt <- grid::gTree(children = grid::gList(ga, gb))
grid::grid.draw(gt)

ga2 <- arrangeGrob(
  arrangeGrob(
    grid::grid.text(label = "Resistance", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    grid::grid.text(label = "", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    con_plot, rgr_plot, nrow = 2, heights = c(0.1,1)),
  ggpubr::as_ggplot(ggpubr::get_legend(plot1(tri_pop, tri_sub, 'Trichomes_sizestd'))),
                    widths = c(1,0.16))
#gc2 <- grid::grid.text(label = "Resistance", hjust = 0, vjust = 1, 
#                x = unit(0.02, "npc"), y = unit(0.98, "npc"), 
#                gp = grid::gpar(fontsize = 15, fontface = "bold"))
gt2 <- grid::gTree(children = grid::gList(ga2, gb))

grid.arrange(gt, gt2, nrow = 2, heights = c(2, 1.1))

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/def-res-through-time-plot.pdf", width = 7.57*1.3, height = 9.2*1.1)
grid.arrange(gt, gt2, nrow = 2, heights = c(2, 1))
dev.off()
png(file = "figures/def-res-through-time-plot.png", width = 7.57*1.3, height = 9.2*1.1,
    units = "in", res = 600)
grid.arrange(gt, gt2, nrow = 2, heights = c(2, 1))
dev.off()

#When finished, remove all the objects used for plotting:
rm(tri_pop, tri_sub, tri_pred, tri_pred2,
   gna_pop, gna_sub, gna_pred, gna_pred2,
   H10_pop, H10_sub, H10_pred, H10_pred2,
   H36_pop, H36_sub, H36_pred, H36_pred2,
   con_pop, con_sub, con_pred, con_pred2,
   rgr_pop, rgr_sub, rgr_pred, rgr_pred2)


kcv_plot <- ggplot() +
    geom_vline(xintercept = -Inf) +
    theme_bw() + theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      axis.title.x = ggtext::element_markdown(size = 12.5),
      axis.title.y = ggtext::element_markdown(size = 12.5),
      axis.text = element_text(size = 12, color = "black"),
      plot.tag.location = "panel",
      plot.tag = element_text(margin = margin(0,0,0,4)),
      plot.margin = margin(10,10,10,10),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, size = 12, color = "black")
    ) +
    labs(x = "Experimental *S. exigua*<br>density (per plant)") +
  scale_y_continuous(limits = c(0,1)) #+
  facet_grid(. ~ Day, labeller = labeller(Day = function(x) paste("Day", x)))

grid.arrange(
  kcv_plot %+%
    geom_point(data = out_sum[which(out_sum$con.notNA >= 7),], 
               aes(x = Herb_add, y = con.KCV),
               size = 1.5, color = "black") +
    # geom_ribbon(data = con_pred, 
    #             aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
    #             fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
    # geom_line(data = con_pred, aes(x = Herb_add, y = median),
    #           color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
    labs(y = "<sup>K</sup>CV of bioassay caterpillar<br>leaf area consumed"),
  
kcv_plot %+%
  geom_point(data = out_sum[which(out_sum$rgr.notNA >= 7),], 
             aes(x = Herb_add, y = rgr.KCV),
             size = 1.5, color = "black") +
  geom_ribbon(data = rgr_pred,
            aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
            fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
  geom_line(data = rgr_pred, aes(x = Herb_add, y = median),
            color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
  labs(y = "<sup>K</sup>CV of bioassay<br>caterpillar RGR"),
# 
# kcv_plot %+%
#   geom_point(data = out_sum, aes(x = Herb_add, y = gna.KCV),
#              size = 1.5, color = "black") +
#   geom_ribbon(data = gna_pred,
#               aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
#               fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
#   geom_line(data = gna_pred, aes(x = Herb_add, y = median),
#             color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
#   labs(y = "<sup>K</sup>CV of gland area"),
# 
# kcv_plot %+%
#   geom_point(data = out_sum, aes(x = Herb_add, y = tri.KCV),
#              size = 1.5, color = "black") +
#   geom_ribbon(data = tri_pred,
#               aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
#               fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
#   geom_line(data = tri_pred, aes(x = Herb_add, y = median),
#             color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
#   labs(y = "<sup>K</sup>CV of trichome density"),
# 
# kcv_plot %+%
#   geom_point(data = out_sum, aes(x = Herb_add, y = H10.KCV),
#              size = 1.5, color = "black") +
#   geom_ribbon(data = H10_pred,
#               aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
#               fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
#   geom_line(data = H10_pred, aes(x = Herb_add, y = median),
#             color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
#   labs(y = "<sup>K</sup>CV of HGQ 10.5"),

kcv_plot %+%
  geom_point(data = out_dam[which(out_dam$Sampling_Percent_damage.notNA >= 7),], 
             aes(x = Herb_add, y = Sampling_Percent_damage.KCV),
             size = 1.5, color = "black") +
  geom_ribbon(data = dam_pred,
              aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
              fill = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], alpha = 0.15) +
  geom_line(data = dam_pred, aes(x = Herb_add, y = median),
            color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[3], linewidth = 1) +
  labs(y = "<sup>K</sup>CV of leaf damage<br>in the field"),
nrow = 1)

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/kcv-plot.pdf", width = 9.72, height = 3)

dev.off()
png(file = "figures/kcv-plot.png", width = 9.72, height = 3,
    units = "in", res = 600)

dev.off()

out_sum[,c('Sampling_Percent_damage.KCV', 'Harvest_Percent_damage.KCV',
           'Plant_RGR.KCV')] <- out_dam[match(out_sum$Cage_pos, out_dam$Plot_ID),
                                        c('Sampling_Percent_damage.KCV', 'Harvest_Percent_damage.KCV',
                                          'Plant_RGR.KCV')]

require(reshape2)
out_sum_melt <- melt(cbind(Cage_pos = out_sum[,'Cage_pos'], out_sum[,grepl("KCV", colnames(out_sum))]),
                     id.vars = 'Cage_pos')
ggplot(data = out_sum_melt, aes(x = variable, y = value)) +
  geom_point()

ggplot(data = out_sum, aes(x = con.KCV, 
                           y = Sampling_Percent_damage.KCV)) +
  geom_point()

rgr_pred <- read.csv("output/kcv-rgr-summary-min5.csv")
rgr_pred <- rgr_pred[grepl("pred", rgr_pred$X),]
#rgr_pred$Timepoint <- as.numeric(with(rgr_pred, substr(X, start = 6, stop = 6)))
rgr_pred$Herb_dens <- as.numeric(with(rgr_pred, substr(X, start = 6, stop = 6)))
#rgr_pred$Day <- c(2, 4, 7, 12, 21)[rgr_pred[,'Timepoint']]
rgr_pred$Herb_add <- c(0, 6, 12, 25)[rgr_pred[,'Herb_dens']]
#rgr_pred <- subset(rgr_pred, Day == 7 | Day == 12)

con_pred <- read.csv("output/kcv-con-summary-min5.csv")
con_pred <- con_pred[grepl("pred", con_pred$X),]
#con_pred$Timepoint <- as.numeric(with(con_pred, substr(X, start = 6, stop = 6)))
con_pred$Herb_dens <- as.numeric(with(con_pred, substr(X, start = 6, stop = 6)))
#con_pred$Day <- c(2, 4, 7, 12, 21)[con_pred[,'Timepoint']]
con_pred$Herb_add <- c(0, 6, 12, 25)[con_pred[,'Herb_dens']]
#con_pred <- subset(con_pred, Day == 7 | Day == 21)

tri_pred <- read.csv("output/kcv-tri-summary-min5.csv")
tri_pred <- tri_pred[grepl("pred", tri_pred$X),]
tri_pred$Herb_dens <- as.numeric(with(tri_pred, substr(X, start = 6, stop = 6)))
tri_pred$Herb_add <- c(0, 6, 12, 25)[tri_pred[,'Herb_dens']]

H10_pred <- read.csv("output/kcv-H10-summary-min5.csv")
H10_pred <- H10_pred[grepl("pred", H10_pred$X),]
H10_pred$Herb_dens <- as.numeric(with(H10_pred, substr(X, start = 6, stop = 6)))
H10_pred$Herb_add <- c(0, 6, 12, 25)[H10_pred[,'Herb_dens']]

gna_pred <- read.csv("output/kcv-gland-area-summary-min5.csv")
gna_pred <- gna_pred[grepl("pred", gna_pred$X),]
gna_pred$Herb_dens <- as.numeric(with(gna_pred, substr(X, start = 6, stop = 6)))
gna_pred$Herb_add <- c(0, 6, 12, 25)[gna_pred[,'Herb_dens']]

dam_pred <- read.csv("output/kcv-damage-summary-min5.csv")
dam_pred <- dam_pred[grepl("pred", dam_pred$X),]
#dam_pred$Timepoint <- as.numeric(with(dam_pred, substr(X, start = 6, stop = 6)))
dam_pred$Herb_dens <- as.numeric(with(dam_pred, substr(X, start = 6, stop = 6)))
#dam_pred$Day <- c(2, 4, 7, 12, 21)[dam_pred[,'Timepoint']]
dam_pred$Herb_add <- c(0, 6, 12, 25)[dam_pred[,'Herb_dens']]
#dam_pred <- subset(dam_pred, Day == 7 | Day == 12)

## =============================================================================
## Strength of induced responses through time (grouped by resistance, physical
##  defenses, chemical defenses)
## =============================================================================

# Read in summary of monitored variables from scaled individual response models.
rgr_sum_sc <- read.csv("output/scaled-rgr-summary.csv")
con_sum_sc <- read.csv("output/scaled-con-summary.csv")
tri_sum_sc <- read.csv("output/scaled-tri-summary.csv")
H10_sum_sc <- read.csv("output/scaled-H10-summary.csv")
H36_sum_sc <- read.csv("output/scaled-H36-summary.csv")
gnd_sum_sc <- read.csv("output/scaled-gland-dens-summary.csv")
gna_sum_sc <- read.csv("output/scaled-gland-area-summary.csv")

# Combine data, adding column for response ID.
trait_sum <- rbind(cbind(rgr_sum_sc, Response = "rgr"),
                   cbind(con_sum_sc, Response = "con"),
                   cbind(tri_sum_sc, Response = "tri"),
                   cbind(H10_sum_sc, Response = "H10"),
                   cbind(H36_sum_sc, Response = "H36"),
                   cbind(gnd_sum_sc, Response = "gnd"),
                   cbind(gna_sum_sc, Response = "gna")
)

trait_pred <- trait_sum[which(sapply(trait_sum$X, FUN = function(x) {grepl("pred", x)})),]
trait_pred$Timepoint <- with(trait_pred, as.numeric(substr(X, start = 6, stop = 6)))
trait_pred$Herb_add <- with(trait_pred, as.numeric(substr(X, start = 8, stop = 8)))
trait_pred$Days <- c(2, 4, 7, 12, 21)[trait_pred[,'Timepoint']]

## Function for plotting size standardized plant level responses with population
##  means (+/- 89% CI). Response ~ days (x axis) x herbivore density (dodged points).
ind_plot <- function(title_name = "", label_vec, color_vec) {
  ggplot(data = subset(trait_pred, Herb_add == 3 & (Response == "rgr" | Response == "con")),
       aes(x = Days, ymin = 0 - CI_high_89, ymax = 0 - CI_low_89, group = Response)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4) +
  geom_ribbon(aes(fill = Response), alpha = 0.1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "inside",
    legend.position.inside = c(0,0),
    legend.justification = c(0,0),
    legend.background = element_blank(),
    axis.title = ggtext::element_markdown(size = 12),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "Days following caterpillar introduction",
       y = "", color = "", fill = "", shape = "", linetype = "",
       title = title_name) +
  scale_x_continuous(breaks = c(0,2,4,7,12,21), limits = c(0, 22), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1.7,2)) +
  scale_color_manual(values = color_vec, labels = label_vec) +
  scale_fill_manual(values = color_vec, labels = label_vec, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed"), labels = label_vec) +
  scale_shape_manual(values = c(16, 17), labels = label_vec)
}

grid.arrange(
ind_plot(title_name = expression(paste(bold(A.), " Induced resistance")),
         label_vec = c("Leaf area consumed", "Bioassay caterpillar RGR"),
         color_vec = RColorBrewer::brewer.pal(n = 8, name = "Accent")[1:2]) %+%
  geom_point(aes(x = Days, y = 0 - mean, color = Response, shape = Response), size = 2.3) +
  geom_line(aes(x = Days, y = 0 - mean, color = Response, linetype = Response), size = 0.8),


ind_plot(title_name = expression(paste(bold(B.), " Induced chemical defenses")),
         label_vec = c("Hemigossypolone 10.5", "Hemigossypolone 3.6"),
         color_vec = RColorBrewer::brewer.pal(n = 8, name = "Accent")[c(3,5)]) %+%
  subset(trait_pred, Herb_add == 3 & (Response == "H10" | Response == "H36")) +
  aes(ymin = CI_low_89, ymax = CI_high_89, group = Response) +
  geom_point(aes(x = Days, y = mean, color = Response, shape = Response), size = 2.3) +
  geom_line(aes(x = Days, y = mean, color = Response, linetype = Response), size = 0.8),

ind_plot(title_name = expression(paste(bold(C.), " Induced physical defenses")),
         label_vec = c("Gland prop. leaf area", "Trichome density"),
         color_vec = RColorBrewer::brewer.pal(n = 8, name = "Accent")[6:7]) %+%
  subset(trait_pred, Herb_add == 3 & (Response == "gna" | Response == "tri")) +
  aes(ymin = CI_low_89, ymax = CI_high_89, group = Response) +
  geom_point(aes(x = Days, y = mean, color = Response, shape = Response), size = 2.3) +
  geom_line(aes(x = Days, y = mean, color = Response, linetype = Response), size = 0.8),

ncol = 1)


dev.off()
# Save (use dev.off() until 'null device 1' comes up)
pdf(file = "figures/induced-plots-1.pdf", width = 4.3, height = 12.25)
#run plot code
dev.off()

dev.off()
# Save (use dev.off() until 'null device 1' comes up)
pdf(file = "figures/induced-plots-2.pdf", width = 7.58, height = 6.35)
#run plot code
dev.off()



################################################################################

# Starting to look at relationship between defenses and resistance at plant and 
#  population levels.
# Work on standardizing each response to respective means of other predictors.

plant_sum <- read.csv("output/plant-con-defen-summary.csv")
full_data <- read.csv("cleaned-data/full-data.csv")

# Get partial residuals 
full_data$Con_defen_resid <- with(full_data, log(Mean_area_consumed) -
                                    (plant_sum[Timepoint, 'mean'] +
                                       plant_sum[which(plant_sum$X == 'slope_tri'), 'mean']*sqrt(Trichomes_std_cm2) +
                                       plant_sum[which(plant_sum$X == 'slope_gnd'), 'mean']*sqrt(Gland_prop) +
                                       plant_sum[which(plant_sum$X == 'slope_H10'), 'mean']*log(H_10.5_relabun) +
                                       plant_sum[which(plant_sum$X == 'slope_H36'), 'mean']*log(H_3.6_relabun) +
                                       plant_sum[which(plant_sum$X == 'slope_resp'), 'mean']*Midvein_length_sc)
                                  )
hist(full_data$Con_defen_resid)
mean(full_data$Con_defen_resid, na.rm = T)

full_data$Trichomes_partial <- with(full_data, Con_defen_resid + 
                                      plant_sum[which(plant_sum$X == 'slope_tri'), 'mean']*sqrt(Trichomes_std_cm2)
                                    )

full_data$Gland_prop_partial <- with(full_data, Con_defen_resid + 
                                      plant_sum[which(plant_sum$X == 'slope_gnd'), 'mean']*sqrt(Gland_prop)
                                     )

full_data$H10.5_partial <- with(full_data, Con_defen_resid + 
                                      plant_sum[which(plant_sum$X == 'slope_H10'), 'mean']*log(H_10.5_relabun)
                                )

full_data$H3.6_partial <- with(full_data, Con_defen_resid + 
                                  plant_sum[which(plant_sum$X == 'slope_H36'), 'mean']*log(H_3.6_relabun)
)


grid.arrange(

ggplot(data = full_data[which(full_data$Mean_area_consumed > 0),], 
       aes(x = sqrt(Trichomes_std_cm2), y = Trichomes_partial)) +
  geom_point() +
  geom_abline(slope = plant_sum[which(plant_sum$X == 'slope_tri'), 'mean']) +
  stat_smooth(method = "lm", color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[4]) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "sqrt Trichomes cm<sup>-2</sup>",
       y = "Leaf area consumed<br>partial residuals"),

ggplot(data = full_data[which(full_data$Mean_area_consumed > 0),], 
       aes(x = sqrt(Gland_prop), y = Gland_prop_partial)) +
  geom_point() +
  geom_abline(slope = plant_sum[which(plant_sum$X == 'slope_gnd'), 'mean']) +
  stat_smooth(method = "lm", color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[4]) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "sqrt Prop. gland area",
       y = "Leaf area consumed<br>partial residuals"),

ggplot(data = full_data[which(full_data$Mean_area_consumed > 0),], 
       aes(x = log(H_3.6_relabun), y = H3.6_partial)) +
  geom_point() +
  #  geom_abline(slope = plant_sum[which(plant_sum$X == 'slope_H10'), 'mean'])
  #  stat_smooth(method = "lm", color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[4]) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "ln Hemigossypolone 3.6",
       y = "Leaf area consumed<br>partial residuals"),

ggplot(data = full_data[which(full_data$Mean_area_consumed > 0),], 
       aes(x = log(H_10.5_relabun), y = H10.5_partial)) +
  geom_point() +
  geom_abline(slope = plant_sum[which(plant_sum$X == 'slope_H10'), 'mean']) +
  stat_smooth(method = "lm", color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[4]) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "ln Hemigossypolone 10.5",
       y = "Leaf area consumed<br>partial residuals") +
  scale_x_continuous(breaks = c(19,20,21)),

nrow = 1 )


dev.off()
pdf(file = "figures/plant-resid-plots.pdf", width = 10.4, height = 3)
#run plot code
dev.off()
png(file = "figures/plant-resid-plots.png", width = 10.4, height = 3,
    units = "in", res = 600)
#run plot code
dev.off()


ggplot(data = full_data[which(full_data$Mean_area_consumed > 0),], 
       aes(x = log(H_3.6_relabun), y = H3.6_partial)) +
  geom_point() +
  #  geom_abline(slope = plant_sum[which(plant_sum$X == 'slope_H10'), 'mean'])
#  stat_smooth(method = "lm", color = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[4]) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "ln Hemigossypolone 3.6",
       y = "Leaf area consumed<br>partial residuals")

dev.off()
pdf(file = "figures/plant-resid-H3.6-plot.pdf", width = 10.4/3, height = 3)
#run plot code
dev.off()
png(file = "figures/plant-resid-H3.6-plot.png", width = 10.4/3, height = 3,
    units = "in", res = 600)
#run plot code
dev.off()



ggplot(data = full_out[-which(full_out$Timepoint == 5),], 
       aes(x = gna, y = rgr, fill = as.factor(Herb_add),
           color = as.factor(Herb_add), shape = as.factor(Herb_add))) +
  geom_errorbar(aes(x = gna, ymin = rgr_CI_low, ymax = rgr_CI_high), alpha = 0.7) +
  geom_errorbarh(aes(xmin = gna_CI_low, xmax = gna_CI_high, y = rgr), alpha = 0.7) +
  geom_point(size = 2.1,
             color = "gray50", stroke = 0.1) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.key.height = unit(0.8, 'cm'),
    legend.key.width = unit(0.8, 'cm'),
    legend.text = element_text(size = 12.5),
    legend.title = element_text(size = 12.5)
  )   +
  #   scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  guides(color = guide_legend(override.aes = list(size = 3.5))) +
  labs(x = "sqrt Prop. gland area (scaled)",
       y = "Bioassay caterpillar RGR (scaled)",
       color = "Experimental<br>*S. exigua*<br>density<br>(per plant)",
       fill = "Experimental<br>*S. exigua*<br>density<br>(per plant)",
       shape = "Experimental<br>*S. exigua*<br>density<br>(per plant)") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_shape_manual(values = c(23, 22, 24, 21))


dev.off()
pdf(file = "figures/gland-area-rgr-plot.pdf", width = 5.27, height = 3.8)
#run plot code
dev.off()
png(file = "figures/gland-area-rgr-plot.png", width = 5.27, height = 3.8,
    units = "in", res = 600)
#run plot code
dev.off()


ggplot(data = full_out[-which(full_out$Timepoint == 5),], 
       aes(x = H10, y = con, fill = as.factor(Herb_add),
           color = as.factor(Herb_add), shape = as.factor(Herb_add))) +
  geom_errorbar(aes(x = H10, ymin = con_CI_low, ymax = con_CI_high), alpha = 0.7) +
  geom_errorbarh(aes(xmin = H10_CI_low, xmax = H10_CI_high, y = con), alpha = 0.7) +
  geom_point(size = 2.1,
             color = "gray50", stroke = 0.1) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = ggtext::element_markdown(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.key.height = unit(0.8, 'cm'),
    legend.key.width = unit(0.8, 'cm'),
    legend.text = element_text(size = 12.5),
    legend.title = element_text(size = 12.5)
  )   +
  #   scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  guides(color = guide_legend(override.aes = list(size = 3.5))) +
  labs(x = "ln Hemigossypolone 10.5 (scaled)",
       y = "ln Leaf area consumed (scaled)",
       color = "Experimental<br>*S. exigua*<br>density<br>(per plant)",
       fill = "Experimental<br>*S. exigua*<br>density<br>(per plant)",
       shape = "Experimental<br>*S. exigua*<br>density<br>(per plant)") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_shape_manual(values = c(23, 22, 24, 21))



dev.off()
pdf(file = "figures/H10-con-plot.pdf", width = 5.27, height = 3.8)
#run plot code
dev.off()
png(file = "figures/H10-con-plot.png", width = 5.27, height = 3.8,
    units = "in", res = 600)
#run plot code
dev.off()


pop_sum <- read.csv("output/pop-con-defen-summary-min5.csv")
pop_sum <- pop_sum[grep("respPop|triPop|gndPop|H10Pop|H36Pop", pop_sum$X),
                   c("X", "mean", "sd", "median", "f", "CI_low_89", "CI_high_89",
                     "CI_low_95", "CI_high_95")]
pop_sum <- pop_sum[-(1:4),]
pop_sum$type <- substr(pop_sum$X, start = 1, stop = 6)

pop_sum_new <- cbind(
  paste0("respPop_", colnames(pop_sum))
)

pop_sum_new <- pop_sum[which(pop_sum$type == "respPo"),2:9]
colnames(pop_sum_new) <- paste0("respPop_", colnames(pop_sum)[2:9])
pop_sum_new[,paste0("triPop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "tri_re"), 2:9]
pop_sum_new[,paste0("gndPop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "gnd_re"), 2:9]
pop_sum_new[,paste0("H10Pop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "H10_re"), 2:9]
pop_sum_new[,paste0("H36Pop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "H36_re"), 2:9]

pop_sum2 <- read.csv("output/pop-con-defen-pop-means-min5.csv")
pop_sum_new[,c("Pop_ID", "Day", "Timepoint", "Herb_add", "Cage_pos")] <-
  pop_sum2[,c("Pop_ID", "Day", "Timepoint", "Herb_add", "Cage_pos")]

pop_fig <- ggplot(data = pop_sum_new,#[-which(pop_sum_new$Timepoint == 5),], 
       aes(x = triPop_mean, y = respPop_mean, color = as.factor(Herb_add))) +
#  geom_errorbar(aes(x = H10Pop_mean, ymin = respPop_CI_low_95, ymax = respPop_CI_high_95),
#                width = 0, linewidth = 0.4) +
#  geom_errorbarh(aes(xmin = H10Pop_CI_low_95, xmax = H10Pop_CI_high_95, y = respPop_mean),
#                 height = 0, linewidth = 0.4) +
  geom_point(size = 2) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = ggtext::element_markdown(size = 12.5),
    axis.title.y = ggtext::element_markdown(size = 12.5),
    axis.text = element_text(size = 12, color = "black"),
    legend.key.height = unit(0.8, 'cm'),
    legend.key.width = unit(0.8, 'cm'),
    legend.text = element_text(size = 12),
    legend.title = ggtext::element_markdown(size = 12, hjust = 0.5),
    plot.tag.location = "panel",
    plot.tag = element_text(margin = margin(0,0,0,4)),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +  
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  labs(x = "Population mean<br>trichome density (scaled)",
       y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
       color = "Experimental<br>*S. exigua*<br>density<br>(per plant)", tag = "A")


g_pt1 <- arrangeGrob(
  g1t,
  arrangeGrob(
pop_fig +
  geom_errorbar(aes(x = triPop_mean, ymin = respPop_CI_low_95, ymax = respPop_CI_high_95),
                width = 0, linewidth = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = triPop_CI_low_95, xmax = triPop_CI_high_95, y = respPop_mean),
                 height = 0, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2) +
  theme(legend.position = "none"),
g1ta, nrow = 2, heights = c(1, 0.2)),

arrangeGrob(
pop_fig %+% aes(x = gndPop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
  geom_errorbar(aes(x = gndPop_mean, ymin = respPop_CI_low_95, ymax = respPop_CI_high_95),
                width = 0, linewidth = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = gndPop_CI_low_95, xmax = gndPop_CI_high_95, y = respPop_mean),
                 height = 0, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2) +
  labs(x = "Population mean<br>prop. gland area (scaled)",
       y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
       color = "Experimental<br>*S. exigua*<br>density<br>(per plant)", tag = "B") +
  theme(legend.position = "none"),
g1ta, nrow = 2, heights = c(1, 0.2)),
nrow = 1,
widths = c(0.2,1,1)
)

g_pt2 <- arrangeGrob(
  g1t,
 arrangeGrob(
  pop_fig %+% aes(x = H36Pop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
    geom_errorbar(aes(x = H36Pop_mean, ymin = respPop_CI_low_95, ymax = respPop_CI_high_95),
                  width = 0, linewidth = 0.4, alpha = 0.5) +
    geom_errorbarh(aes(xmin = H36Pop_CI_low_95, xmax = H36Pop_CI_high_95, y = respPop_mean),
                   height = 0, linewidth = 0.4, alpha = 0.5) +
    geom_point(size = 2) +
    labs(x = "Population mean<br>HGQ 3.6 (scaled)",
         y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
         color = "Experimental<br>*S. exigua*<br>density<br>(per plant)", tag = "C") +
    theme(legend.position = "none"),
  g1ta, nrow = 2, heights = c(1, 0.2)),
  
 arrangeGrob(
  pop_fig %+% aes(x = H10Pop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
    geom_errorbar(aes(x = H10Pop_mean, ymin = respPop_CI_low_95, ymax = respPop_CI_high_95),
                  width = 0, linewidth = 0.4, alpha = 0.5) +
    geom_errorbarh(aes(xmin = H10Pop_CI_low_95, xmax = H10Pop_CI_high_95, y = respPop_mean),
                   height = 0, linewidth = 0.4, alpha = 0.5) +
    geom_point(size = 2) +
    labs(x = "Population mean<br>HGQ 10.5 (scaled)",
         y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
         color = "Experimental<br>*S. exigua*<br>density<br>(per plant)", tag = "D")+
    theme(legend.position = "none"),
  g1ta, nrow = 2, heights = c(1, 0.2)),
  nrow = 1,
  widths = c(0.2,1,1)
)



g_pts <- arrangeGrob(g_pt1, g_pt2, nrow = 2)
grid::grid.newpage()
grid.arrange(g_pts, ggpubr::as_ggplot(ggpubr::get_legend(pop_fig)),
             ncol = 2, widths = c(1, 0.15))

# Save (use dev.off() until 'null device 1' comes up)
dev.off()
pdf(file = "figures/pop-def-res-plot.pdf", width = 8.5, height = 7.8*0.96)
grid.arrange(g_pts, ggpubr::as_ggplot(ggpubr::get_legend(pop_fig)),
             ncol = 2, widths = c(1, 0.15))
dev.off()
png(file = "figures/pop-def-res-plot.png", width = 8.5, height = 7.8*0.96,
    units = "in", res = 600)
grid.arrange(g_pts, ggpubr::as_ggplot(ggpubr::get_legend(pop_fig)),
             ncol = 2, widths = c(1, 0.15))
dev.off()


g1 <- grid::grid.move.to(x = 0.75, y = 0.95)
g2 <- grid::grid.line.to(x = 0.75, y = 0.35, arrow = grid::arrow(),
                   gp = grid::gpar(lwd = 3))
g3 <- grid::grid.text(label = "More\nresistant", hjust = 0.5, vjust = 0.5, 
                x = unit(0.75, "npc"), y = unit(0.25, "npc"), 
                gp = grid::gpar(fontsize = 10, fontface = "bold.italic"))
g1t <- grid::gTree(children = grid::gList(g1, g2, g3), gp = grid::gpar(hjust = 0))
g1a <- grid::grid.move.to(x = 0.25, y = 0.65)
g2a <- grid::grid.line.to(x = 0.8, y = 0.65, arrow = grid::arrow(),
                         gp = grid::gpar(lwd = 3))
g3a <- grid::grid.text(label = "More\ndefended", hjust = 0.5, vjust = 0.5, 
                      x = unit(0.9, "npc"), y = unit(0.65, "npc"), 
                      gp = grid::gpar(fontsize = 10, fontface = "bold.italic"))
g1ta <- grid::gTree(children = grid::gList(g1a, g2a, g3a), gp = grid::gpar(hjust = 0))

full_data[which(full_data$Gland_prop_sizestd > 0.24 &
                  full_data$Gland_prop_sizestd < 0.26 &
                  full_data$Midvein_length > 3.8 &
                  full_data$Midvein_length < 4.2),'Plant_ID']

full_data[which(full_data$Gland_prop_sizestd > 0.31 &
                  full_data$Gland_prop_sizestd < 0.33 &
                  full_data$Midvein_length > 3.8 &
                  full_data$Midvein_length < 4.2),'Plant_ID']

full_data[which(full_data$Gland_prop_sizestd > 0.16 &
                  full_data$Gland_prop_sizestd < 0.19 &
                  full_data$Midvein_length > 3.8 &
                  full_data$Midvein_length < 4.2),'Plant_ID']


hist(full_data$Gland_prop_sizestd)
