## USDA cotton experiment
## Figures
## Code written and annotated by JM (last updated 2026-03-27)

# The code in this script file:
# - Produces Figures 1-4 from the main text.

# Load required packages.
require(ggplot2)
require(gridExtra)
require(ggtext)

## =============================================================================
## Figure 1. Ant attack as a function of experimental S. exigua density.
## =============================================================================

# Read in ant attack data.
ants <- read.csv("cleaned-data/ant-attack.csv")

# Plot
ggplot() +
  geom_jitter(data = ants, aes(x = Herb_add, y = No_plants_attacked/No_plants_obs),
              height = 0, width = 1, size = 1.5, alpha = 0.6) +
  geom_line(aes(x = c(6,12,25), y = c(0.11,0.18,0.4)), size = 1, color = "darkred") +
  geom_point(aes(x = c(6,12,25), y = c(0.11,0.18,0.4)), 
             size = 3, color = "darkred") +
  geom_errorbar(aes(x = c(6,12,25), ymin = c(0.07,0.14,0.33), ymax = c(0.15,0.22,0.47)),
                width = 0.6, color = "darkred") +
  
  theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x = ggtext::element_markdown(size = 12.5),
    axis.title.y = ggtext::element_markdown(size = 12.5),
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "*Spodoptera* plant<sup>-1</sup>", 
       y = "Prop. plants tended by ants") +
  scale_y_continuous(limits = c(0,1)) 

# # Save figure.
dev.off()
pdf(file = "figures/ant-attack-plot-2026-10-08.pdf", width = 4.2, height = 3.3)
#run plot code
dev.off()
png(file = "figures/ant-attack-plot-2026-10-08.png", width = 4.2, height = 3.3,
    units = "in", res = 600)
#run plot code
dev.off()

## =============================================================================
## Figure 2. Plant responses with populations means (across timepoints and herbivore densities)
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
gna_summary[which(gna_summary$X == "slope_resp"), ]

# For response variables with significant effect of leaf size, standardize trait 
#  values to mean leaf size
full_data$rgr_sizestd <- with(full_data, RGR_mean - 
                                rgr_summary[which(rgr_summary$X == "slope_resp"), 'mean']*Initial_leaf_area_sc)
full_data$con_sizestd <- with(full_data, ifelse(Mean_area_consumed == 0, NA,
                                                log(Mean_area_consumed)))
# unscale leaf area consumed (JM 2026-07-27)
full_data$con_sizestd_unsc <- with(full_data, exp(con_sizestd))
full_data$Trichomes_sizestd <- with(full_data, sqrt(Trichomes_std_cm2) - 
                                      tri_summary[which(tri_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$Trichomes_sizestd_unsc <- with(full_data, Trichomes_std_cm2 - 
                                      (tri_summary[which(tri_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)^2)
# unscale gland proportional area (JM 2026-07-26)
full_data$H_10.5_sizestd <- with(full_data, log(H_10.5_relabun))
full_data$H_3.6_sizestd <- with(full_data, log(H_3.6_relabun) - 
                                  H36_summary[which(H36_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$Gland_prop_sizestd <- with(full_data, sqrt(Gland_prop) - 
                                       gna_summary[which(gna_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
# unscale gland proportional area (JM 2026-07-26)
full_data$Gland_prop_sizestd_unsc <- with(full_data, Gland_prop - 
                                       (gna_summary[which(gna_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)^2)


# Read in population means for each response.
rgr_pop <- read.csv("output/full-rgr-pop-means-min5.csv")
con_pop <- read.csv("output/full-con-pop-means-min5.csv")
tri_pop <- read.csv("output/full-tri-pop-means-min5.csv")
gnd_pop <- read.csv("output/full-gland-dens-pop-means-min5.csv")
gna_pop <- read.csv("output/full-gland-area-pop-means-min5.csv")
H10_pop <- read.csv("output/full-H10-pop-means-min5.csv")
H36_pop <- read.csv("output/full-H36-pop-means-min5.csv")

## Function for plotting size standardized plant level responses with population
##  means (+/- 95% CI). Response ~ days (x axis) x herbivore density (dodged points).
plot1 <- function(pop_df, model_df, pop_response, model_response) {
  ggplot(data = pop_df, aes(pop_df[,'Day'], y = pop_df[,pop_response], 
                            group = as.factor(Herb_add),
                            fill = as.factor(Herb_add),
                            shape = as.factor(Herb_add))) +
    geom_point(data = model_df, aes(x = model_df[,'Day'], y = model_df[,model_response],
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
    labs(x = "Days since *Spodoptera* introduction",
         y = "sqrt Trichomes cm<sup>-2</sup>", 
         fill = "*Spodoptera*<br>plant<sup>-1</sup>", 
         shape = "*Spodoptera*<br>plant<sup>-1</sup>",
         color = "*Spodoptera*<br>plant<sup>-1</sup>") +
    guides(color = guide_legend(override.aes = list(size = 3.5))) +
    scale_x_continuous(breaks = c(0,2,4,7,12,21), limits = c(0, 22), expand = c(0,0)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_shape_manual(values = c(23, 22, 24, 21))
}

cols <- c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 'Day')

# ## Trichome density as response.
# # Set up dataframe.
# tri_sub <- full_data[,c(cols, 'Trichomes_sizestd', 'Midvein_length_sc')]
# tri_sub <- na.omit(tri_sub)
# # Include only populations with at least 5 plant observations.
# tri_sub <- tri_sub[which(tri_sub$Pop_ID %in% names(which(table(tri_sub$Pop_ID) >= 5))),]
# 
# # Set up dataframe with herbivore density x timepoint predictions from 'summary'
# tri_pred2 <- cbind(
#   Day = rep(c(2,4,7,12,21),4),
#   Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
#   pred = c(tri_summary[1:5, 'mean'], tri_summary[1:5, 'mean'] + 
#              tri_summary[grepl('pred', tri_summary$X), 'mean'])
# )

### THIS IS THE UNSCALED VERSION
tri_sub <- full_data[,c(cols, 'Trichomes_sizestd_unsc', 'Midvein_length_sc')]
tri_sub <- na.omit(tri_sub)
# Include only populations with at least 5 plant observations.
tri_sub <- tri_sub[which(tri_sub$Pop_ID %in% names(which(table(tri_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
tri_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(tri_summary[1:5, 'mean'], tri_summary[1:5, 'mean'] + 
             tri_summary[grepl('pred', tri_summary$X), 'mean'])
)

tri_pop$mean_unsc <- with(tri_pop, mean^2)
tri_pop$CI_low_95 <- with(tri_pop, CI_low_95^2)
tri_pop$CI_high_95 <- with(tri_pop, CI_high_95^2)

tri_plot <- plot1(pop_df = tri_pop, model_df = tri_sub, 
      pop_response = 'mean_unsc', model_response = 'Trichomes_sizestd_unsc') +
  theme(legend.position = "none")+
  labs(y = "Trichomes cm<sup>-2</sup>", tag = "A") +
  annotate("text", x = c(21), y = c(Inf), vjust = 1,
           label = c("+")) +
  geom_line(data = tri_pred2, aes(x = Day, y = pred^2, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none") 


## Hemigossypolone 10.5 relative abundance as response.
# Set up dataframe.
H10_sub <- full_data[,c(cols, 'H_10.5_sizestd', 'Midvein_length_sc')]
H10_sub <- na.omit(H10_sub)
# Include only populations with at least 5 plant observations.
H10_sub <- H10_sub[which(H10_sub$Pop_ID %in% names(which(table(H10_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
H10_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(H10_summary[1:5, 'mean'], H10_summary[1:5, 'mean'] + 
             H10_summary[grepl('pred', H10_summary$X), 'mean'])
)

## Hemigossypolone 3.6 relative abundance as response.
# Set up dataframe
H36_sub <- full_data[,c(cols, 'H_3.6_sizestd', 'Midvein_length_sc')]
H36_sub <- na.omit(H36_sub)
# Include only populations with at least 5 plant observations.
H36_sub <- H36_sub[which(H36_sub$Pop_ID %in% names(which(table(H36_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
H36_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(H36_summary[1:5, 'mean'], H36_summary[1:5, 'mean'] + 
             H36_summary[grepl('pred', H36_summary$X), 'mean'])
)

# ## Gland proportional leaf area as response.
# # Set up dataframe.
# gna_sub <- full_data[,c(cols, 'Gland_prop_sizestd', 'Midvein_length_sc')]
# gna_sub <- na.omit(gna_sub)
# # Include only populations with at least 5 plant observations.
# gna_sub <- gna_sub[which(gna_sub$Pop_ID %in% names(which(table(gna_sub$Pop_ID) >= 5))),]
# 
# # Set up dataframe with herbivore density x timepoint predictions from 'summary'
# gna_pred2 <- cbind(
#   Day = rep(c(2,4,7,12,21),4),
#   Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
#   pred = c(gna_summary[1:5, 'mean'], gna_summary[1:5, 'mean'] + 
#              gna_summary[grepl('pred', gna_summary$X), 'mean'])
# )


#### THIS IS THE UNSCALED VERSION
## Gland proportional leaf area as response.
# Set up dataframe.
gna_sub <- full_data[,c(cols, 'Gland_prop_sizestd_unsc', 'Midvein_length_sc')]
gna_sub <- na.omit(gna_sub)
# Include only populations with at least 5 plant observations.
gna_sub <- gna_sub[which(gna_sub$Pop_ID %in% names(which(table(gna_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
gna_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(gna_summary[1:5, 'mean'], gna_summary[1:5, 'mean'] + 
             gna_summary[grepl('pred', gna_summary$X), 'mean'])
)

gna_pop$mean_unsc <- with(gna_pop, mean^2)
gna_pop$CI_low_95 <- with(gna_pop, CI_low_95^2)
gna_pop$CI_high_95 <- with(gna_pop, CI_high_95^2)

gna_plot <- plot1(pop_df = gna_pop, model_df = gna_sub, 
      pop_response = 'mean_unsc', model_response = 'Gland_prop_sizestd_unsc') +
  theme(legend.position = "none")+
  labs(y = "Pigment gland prop. leaf area", tag = "B") +
  annotate("text", x = c(7,21), y = c(Inf, Inf), vjust = 1,
           label = c("*", "*")) +
  geom_line(data = gna_pred2, aes(x = Day, y = pred^2, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none") +
  scale_y_continuous(limits = c(0, 0.2))


## Bioassay leaf area consumed as response.
# Set up dataframe.
con_sub <- full_data[,c(cols, 'con_sizestd', 'Initial_leaf_area_sc')]
con_sub <- na.omit(con_sub)
# Include only populations with at least 5 plant observations.
con_sub <- con_sub[which(con_sub$Pop_ID %in% names(which(table(con_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
con_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(con_summary[1:5, 'mean'], con_summary[1:5, 'mean'] + 
             con_summary[grepl('pred', con_summary$X), 'mean'])
)

# ######### THIS IS THE UNSCALED VERSION
# ## Bioassay leaf area consumed as response.
# # Set up dataframe.
# con_sub <- full_data[,c(cols, 'con_sizestd_unsc', 'Initial_leaf_area_sc')]
# con_sub <- na.omit(con_sub)
# # Include only populations with at least 5 plant observations.
# con_sub <- con_sub[which(con_sub$Pop_ID %in% names(which(table(con_sub$Pop_ID) >= 5))),]
# 
# # Set up dataframe with herbivore density x timepoint predictions from 'summary'
# con_pred2 <- cbind(
#   Day = rep(c(2,4,7,12,21),4),
#   Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
#   pred = c(con_summary[1:5, 'mean'], con_summary[1:5, 'mean'] + 
#              con_summary[grepl('pred', con_summary$X), 'mean'])
# )
# 
# con_pop$mean_unsc <- with(con_pop, exp(mean))
# con_pop$CI_low_95 <- with(con_pop, exp(CI_low_95))
# con_pop$CI_high_95 <- with(con_pop, exp(CI_high_95))
# 
# 
#   plot1(pop_df = con_pop, model_df = con_sub,
#         pop_response = 'mean_unsc', model_response = 'con_sizestd_unsc') +
#   labs(y = "Bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)") +
#   annotate("text", x = c(2,4,7,12), y = c(Inf, Inf, Inf, Inf), vjust = 1,
#            label = c("*", "*", "*", "*")) +
#   geom_line(data = con_pred2, aes(x = Day, y = exp(pred), color = as.factor(Herb_add),
#                                   linetype = as.factor(Herb_add)),
#             position = position_dodge(width = 1.5), linewidth = 0.3) +
#   scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"),
#                         guide = "none")
# 
#   # Save (use dev.off() until 'null device 1' comes up)
#   dev.off()
#   pdf(file = "figures/con-through-time-unsc-2026-10-08.pdf", width = 6.5, height = 3.5)
# 
#   dev.off()
#   png(file = "figures/con-through-time-unsc-2026-10-08.png", width = 6.5, height = 3.5,
#       units = "in", res = 600)
# 
#   dev.off()

## Bioassay caterpillar RGR as response.
# Set up dataframe.
rgr_sub <- full_data[,c(cols, 'rgr_sizestd', 'Initial_leaf_area_sc')]
rgr_sub <- na.omit(rgr_sub)
# Include only populations with at least 5 plant observations.
rgr_sub <- rgr_sub[which(rgr_sub$Pop_ID %in% names(which(table(rgr_sub$Pop_ID) >= 5))),]

# Set up dataframe with herbivore density x timepoint predictions from 'summary'
rgr_pred2 <- cbind(
  Day = rep(c(2,4,7,12,21),4),
  Herb_add = c(rep(0,5), rep(6,5), rep(12,5), rep(25,5)),
  pred = c(rgr_summary[1:5, 'mean'], rgr_summary[1:5, 'mean'] + 
             rgr_summary[grepl('pred', rgr_summary$X), 'mean'])
)

# Plot panels for each defense or measure of resistance.
# tri_plot <-
#   plot1(pop_df = tri_pop, model_df = tri_sub, 
#         pop_response = 'mean', model_response = 'Trichomes_sizestd') +
#   theme(legend.position = "none") +
#   labs(y = "sqrt Trichomes cm<sup>-2</sup>", tag = "A") +
#   annotate("text", x = c(21), y = c(Inf), vjust = 1, label = c("+")) +
#   geom_line(data = tri_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
#                                   linetype = as.factor(Herb_add)),
#             position = position_dodge(width = 1.5), linewidth = 0.3) +
#   scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
#                         guide = "none")
  
H10_plot <- 
  plot1(pop_df = H10_pop, model_df = H10_sub, 
        pop_response = 'mean', model_response = 'H_10.5_sizestd') +
  theme(legend.position = "none") +
  labs(y = "ln HGQ 10.5", tag = "D") +
  annotate("text", x = c(4,7,12), y = c(Inf, Inf, Inf), vjust = 1, 
           label = c("*", "*", "*")) +
  geom_line(data = H10_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none")
  
H36_plot <- 
  plot1(pop_df = H36_pop, model_df = H36_sub, 
        pop_response = 'mean', model_response = 'H_3.6_sizestd') +
  theme(legend.position = "none") +
  labs(y = "ln HGQ 3.6", tag = "C")  +
  annotate("text", x = c(4,12), y = c(Inf, Inf), vjust = 1, label = c("*", "+")) +
  geom_line(data = H36_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none")

# gna_plot <-
#   plot1(pop_df = gna_pop, model_df = gna_sub, 
#         pop_response = 'mean', model_response = 'Gland_prop_sizestd') +
#   theme(legend.position = "none") +
#   labs(y = "sqrt Prop. gland area", tag = "B") +
#   annotate("text", x = c(7,21), y = c(Inf, Inf), vjust = 1, label = c("*", "*")) +
#   geom_line(data = gna_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
#                                   linetype = as.factor(Herb_add)),
#             position = position_dodge(width = 1.5), linewidth = 0.3) +
#   scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
#                         guide = "none")
  
con_plot <-
  plot1(pop_df = con_pop, model_df = con_sub, 
        pop_response = 'mean', model_response = 'con_sizestd') +
  theme(legend.position = "none")+
  labs(y = "ln Bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)", tag = "E") +
  annotate("text", x = c(2,4,7,12), y = c(Inf, Inf, Inf, Inf), vjust = 1,
           label = c("*", "*", "*", "*")) +
  geom_line(data = con_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none")
  
rgr_plot <-
  plot1(pop_df = rgr_pop, model_df = rgr_sub, 
        pop_response = 'mean', model_response = 'rgr_sizestd') +
  theme(legend.position = "none")+
  labs(y = "Bioassay caterpillar<br>RGR (g g<sup>-1</sup> d<sup>-1</sup>)", tag = "F") +
  annotate("text", x = c(4,7,21), y = c(Inf, Inf, Inf), vjust = 1,
           label = c("*", "*", "+")) +
  geom_line(data = rgr_pred2, aes(x = Day, y = pred, color = as.factor(Herb_add),
                                  linetype = as.factor(Herb_add)),
            position = position_dodge(width = 1.5), linewidth = 0.3) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), 
                        guide = "none")
  
# Get legend for use in arranged grobs.  
ggpubr::as_ggplot(ggpubr::get_legend(plot1(pop_df = tri_pop, model_df = tri_sub, 
                                           pop_response = 'mean',
                                           model_response = 'Trichomes_sizestd_unsc'
                                           )))  

# Make "Defenses" box: "Defenses" label, panels, legend (ga), and outline (gb).
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
  ggpubr::as_ggplot(ggpubr::get_legend(plot1(tri_pop, tri_sub, 'mean', 'Trichomes_sizestd_unsc'))),
                    ncol = 2, widths = c(1,0.16) )
gb <- grid::grid.rect(width = 0.985, height = 0.985, 
                gp = grid::gpar(lwd = 1.5, col = "black", fill = NA))
gt <- grid::gTree(children = grid::gList(ga, gb))
grid::grid.draw(gt)

# Make "Resistance" box: "Resistance" label, panels, legend (ga2), and outline (gb2).
ga2 <- arrangeGrob(
  arrangeGrob(
    grid::grid.text(label = "Resistance", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    grid::grid.text(label = "", hjust = 0, vjust = 1, 
                    x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                    gp = grid::gpar(fontsize = 15, fontface = "bold")),
    con_plot, rgr_plot, nrow = 2, heights = c(0.1,1)),
  ggpubr::as_ggplot(ggpubr::get_legend(plot1(tri_pop, tri_sub, 'mean', 'Trichomes_sizestd_unsc'))),
                    widths = c(1,0.16))
gt2 <- grid::gTree(children = grid::gList(ga2, gb))

# Put them together.
p <- grid.arrange(gt, gt2, nrow = 2, heights = c(2, 1.1))

# Save
ggsave("figures/def-res-through-time-plot-2026-11-02.pdf", plot = p,
       width = 7.57*1.3, height = 9.2*1.1, device = cairo_pdf)
ggsave("figures/def-res-through-time-plot-2026-11-02.png", plot = p,
       width = 7.57*1.3, height = 9.2*1.1, dpi = 600)


#When finished, remove all the objects used for plotting:
rm(p, tri_pop, tri_sub, tri_summary, tri_pred2,
   gna_pop, gna_sub, gna_summary, gna_pred2,
   H10_pop, H10_sub, H10_summary, H10_pred2,
   H36_pop, H36_sub, H36_summary, H36_pred2,
   con_pop, con_sub, con_summary, con_pred2,
   rgr_pop, rgr_sub, rgr_summary, rgr_pred2)

## =============================================================================
## Figure 3. KCV of resistance and damage as a function of experimental S. exigua density.
## =============================================================================

# Read in KCV estimates.
out_sum <- read.csv('output/KCV-estimates.csv')

# Plotting function
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
    labs(x = "*Spodoptera* plant<sup>-1</sup>") +
  scale_y_continuous(limits = c(0,1.07), breaks = seq(0,1,by = 0.25)) 

# For bioassay caterpillar RGR, bioassay leaf area consumed, and leaf damage at
#  sampling, use model summaries to get predictions for KCV of each response at 
#  each herbivore density.
rgr_pred <- read.csv("output/kcv-rgr-summary-min5.csv")
rgr_pred <- rgr_pred[grepl("pred", rgr_pred$X),]
rgr_pred$Herb_dens <- as.numeric(with(rgr_pred, substr(X, start = 6, stop = 6)))
rgr_pred$Herb_add <- c(0, 6, 12, 25)[rgr_pred[,'Herb_dens']]

con_pred <- read.csv("output/kcv-con-summary-min5.csv")
con_pred <- con_pred[grepl("pred", con_pred$X),]
con_pred$Herb_dens <- as.numeric(with(con_pred, substr(X, start = 6, stop = 6)))
con_pred$Herb_add <- c(0, 6, 12, 25)[con_pred[,'Herb_dens']]
  
dam_pred <- read.csv("output/kcv-damage-summary-min5.csv")
dam_pred <- dam_pred[grepl("pred", dam_pred$X),]
dam_pred$Herb_dens <- as.numeric(with(dam_pred, substr(X, start = 6, stop = 6)))
dam_pred$Herb_add <- c(0, 6, 12, 25)[dam_pred[,'Herb_dens']]

# Use plotting function to plot KCV for bioassay leaf area consumed, bioassay
#  caterpillar RGR, and leaf damage at sampling, with model predictions.
p <- grid.arrange(
  kcv_plot %+%
    geom_point(data = out_sum[which(out_sum$con.notNA >= 7),], 
               aes(x = Herb_add, y = con.KCV, 
                   color = as.factor(Herb_add), fill = as.factor(Herb_add),
                   shape = as.factor(Herb_add)),
               size = 1.5) +
    theme(legend.position = "none") +
    labs(y = "<sup>K</sup>CV of bioassay caterpillar<br>leaf area consumed", tag = "A")  +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
    scale_shape_manual(values = c(23, 22, 24, 21)),
  
kcv_plot %+%
  geom_point(data = out_sum[which(out_sum$rgr.notNA >= 7),], 
             aes(x = Herb_add, y = rgr.KCV, 
                 color = as.factor(Herb_add), fill = as.factor(Herb_add),
                 shape = as.factor(Herb_add)),
             size = 1.5) +
  geom_ribbon(data = rgr_pred,
            aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
            fill = "black", alpha = 0.15) +
  geom_line(data = rgr_pred, aes(x = Herb_add, y = median),
            color = "black") +
  theme(legend.position = "none") +
  labs(y = "<sup>K</sup>CV of bioassay<br>caterpillar RGR", tag = "B") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_shape_manual(values = c(23, 22, 24, 21)),

kcv_plot %+%
  geom_point(data = out_sum[which(out_sum$samp_dam.notNA >= 7),], 
             aes(x = Herb_add, y = samp_dam.KCV, 
                 color = as.factor(Herb_add), fill = as.factor(Herb_add),
                 shape = as.factor(Herb_add)),
             size = 1.5) +
  geom_ribbon(data = dam_pred,
              aes(x = Herb_add, ymin = CI_low_95, ymax = CI_high_95),
              fill = "black", alpha = 0.15) +
  geom_line(data = dam_pred, aes(x = Herb_add, y = median),
            color = "black") +
  theme(legend.position = "none") +
  labs(y = "<sup>K</sup>CV of leaf damage<br>in the field", tag = "C") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, name = "GnBu")[2:5]) +
  scale_shape_manual(values = c(23, 22, 24, 21)),
nrow = 1)

# Save
ggsave("figures/kcv-plot.pdf", plot = p,
       width = 9.72, height = 3, device = cairo_pdf)
ggsave("figures/kcv-plot.png", plot = p,
       width = 9.72, height = 3, dpi = 600)

rm(p)

## =============================================================================
## Figure 4. Relationship between population mean defenses and resistance.
## =============================================================================

# Read in model summary from hierarchical model.
pop_sum <- read.csv("output/pop-con-defen-summary-min5.csv")
# Get model estimates for experimental population means -- this will include
#  pop means for all defenses and resistance in long format.
pop_sum <- pop_sum[grep("respPop|triPop|gndPop|H10Pop|H36Pop", pop_sum$X),
                   c("X", "mean", "sd", "median", "f", "CI_low_89", "CI_high_89",
                     "CI_low_95", "CI_high_95")]
pop_sum <- pop_sum[-(1:4),]
# Add column for type of population mean: 'respPo' will be pop mean resistance,
#  'XXX_re' will be pop mean value for defense XXX.
pop_sum$type <- substr(pop_sum$X, start = 1, stop = 6)

# Convert to wide format with one row per experimental population.
# Create new dataframe with one row for each experimental population.
pop_sum_new <- pop_sum[which(pop_sum$type == "respPo"),2:9]
# Rename columns to indicate that these values are for pop mean resistance.
colnames(pop_sum_new) <- paste0("respPop_", colnames(pop_sum)[2:9])
# Add new columns for pop means for each defense.
pop_sum_new[,paste0("triPop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "tri_re"), 2:9]
pop_sum_new[,paste0("gndPop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "gnd_re"), 2:9]
pop_sum_new[,paste0("H10Pop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "H10_re"), 2:9]
pop_sum_new[,paste0("H36Pop_", colnames(pop_sum)[2:9])] <- 
  pop_sum[which(pop_sum$type == "H36_re"), 2:9]

# Read in model output on population means, which has metadata on the actual
#  population ID, herbivore density treatment, sampling timepoint, etc, for
#  each experimental population.
pop_sum2 <- read.csv("output/pop-con-defen-pop-means-min5.csv")
# Add these metadata to the wide format dataframe.
pop_sum_new[,c("Pop_ID", "Day", "Timepoint", "Herb_add", "Cage_pos")] <-
  pop_sum2[,c("Pop_ID", "Day", "Timepoint", "Herb_add", "Cage_pos")]

# Plotting function
pop_fig <- ggplot(data = pop_sum_new,#[-which(pop_sum_new$Timepoint == 5),], 
       aes(x = triPop_mean, y = respPop_mean, color = as.factor(Herb_add))) +
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
  labs(x = "Population&nbsp;mean<br>trichome&nbsp;density&nbsp;(scaled)",
       y = "Population&nbsp;mean<br>ln&nbsp;bioassay&nbsp;caterpillar<br>leaf&nbsp;area&nbsp;consumed (cm<sup>2</sup>)",
       color = "*Spodoptera*<br>plant<sup>-1</sup>", tag = "A")

# Make plotting objects for axis arrows and labels.
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

# Arrange plotting objects with panels for trichome density and gland prop. area.
g_pt1 <- arrangeGrob(
  g1t,
  arrangeGrob(
    pop_fig +
      geom_errorbar(aes(x = triPop_mean, ymin = respPop_CI_low_95, 
                        ymax = respPop_CI_high_95),
                    width = 0, linewidth = 0.4, alpha = 0.5) +
      geom_errorbarh(aes(xmin = triPop_CI_low_95, xmax = triPop_CI_high_95, 
                         y = respPop_mean),
                     height = 0, linewidth = 0.4, alpha = 0.5) +
      geom_point(size = 2) +
      theme(legend.position = "none"),
    g1ta, nrow = 2, heights = c(1, 0.2)),
  
  arrangeGrob(
    pop_fig %+% aes(x = gndPop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
      geom_errorbar(aes(x = gndPop_mean, ymin = respPop_CI_low_95, 
                        ymax = respPop_CI_high_95),
                    width = 0, linewidth = 0.4, alpha = 0.5) +
      geom_errorbarh(aes(xmin = gndPop_CI_low_95, xmax = gndPop_CI_high_95, 
                         y = respPop_mean),
                     height = 0, linewidth = 0.4, alpha = 0.5) +
      geom_point(size = 2) +
      labs(x = "Population mean<br>prop. gland area (scaled)",
           y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
           color = "*Spodoptera*<br>plant<sup>-1</sup>", tag = "B") +
      theme(legend.position = "none"),
    g1ta, nrow = 2, heights = c(1, 0.2)),
  
  nrow = 1, widths = c(0.2,1,1)
)

# Arrange plotting objects with panels for HGQ 3.6 and 10.5.
g_pt2 <- arrangeGrob(
  g1t,
  arrangeGrob(
    pop_fig %+% aes(x = H36Pop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
      geom_errorbar(aes(x = H36Pop_mean, ymin = respPop_CI_low_95, 
                        ymax = respPop_CI_high_95),
                    width = 0, linewidth = 0.4, alpha = 0.5) +
      geom_errorbarh(aes(xmin = H36Pop_CI_low_95, xmax = H36Pop_CI_high_95, 
                         y = respPop_mean),
                     height = 0, linewidth = 0.4, alpha = 0.5) +
      geom_point(size = 2) +
      labs(x = "Population mean<br>HGQ 3.6 (scaled)",
           y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
           color = "*Spodoptera*<br>plant<sup>-1</sup>", tag = "C") +
      theme(legend.position = "none"),
    g1ta, nrow = 2, heights = c(1, 0.2)),
  
  arrangeGrob(
    pop_fig %+% aes(x = H10Pop_mean, y = respPop_mean, color = as.factor(Herb_add)) +
      geom_errorbar(aes(x = H10Pop_mean, ymin = respPop_CI_low_95, 
                        ymax = respPop_CI_high_95),
                    width = 0, linewidth = 0.4, alpha = 0.5) +
      geom_errorbarh(aes(xmin = H10Pop_CI_low_95, xmax = H10Pop_CI_high_95, 
                         y = respPop_mean),
                     height = 0, linewidth = 0.4, alpha = 0.5) +
      geom_point(size = 2) +
      labs(x = "Population mean<br>HGQ 10.5 (scaled)",
           y = "Population mean<br>ln bioassay caterpillar<br>leaf area consumed (cm<sup>2</sup>)",
           color = "*Spodoptera*<br>plant<sup>-1</sup>", tag = "D")+
      theme(legend.position = "none"),
    g1ta, nrow = 2, heights = c(1, 0.2)),
  
  nrow = 1, widths = c(0.2,1,1)
)

# Put them together.
g_pts <- arrangeGrob(g_pt1, g_pt2, nrow = 2)
grid::grid.newpage()
p <- grid.arrange(g_pts, ggpubr::as_ggplot(ggpubr::get_legend(pop_fig)),
             ncol = 2, widths = c(1, 0.15))

# Save
ggsave("figures/pop-def-res-plot 2026-11-02.pdf", plot = p,
       width = 8.5, height = 7.8*0.96, device = cairo_pdf)
ggsave("figures/pop-def-res-plot 2026-11-02.png", plot = p,
       width = 8.5, height = 7.8*0.96, dpi = 600)
