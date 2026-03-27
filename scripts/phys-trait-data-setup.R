## USDA cotton experiment
## Physical trait data setup
## Code written and annotated by JM (last updated 2024-10-16)

# The code in this script file:
# - standardizes trichome and terpenoid gland counts by observer;
# - converts trichome and terpenoid gland counts to densities per sq cm and 
#   averages these densities by plant;
# - combines average trichome and terpenoid densities with terpenoid gland area
#   and midvein length of sampled leaf;
# - saves the cleaned dataframe as .csv.

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("./USDA-cotton-exp")

# Load required packages.
require(ggplot2)
require(doBy)

## First, use subset of leaf images measured by more than one observer to 
##  standardize trichome and gland counts across observers.

# Read in standardization dataset
phys_traits_obs <- read.csv("raw-data/phys-traits-by-observer.csv")

# View relationship between gland counts by JM and gland counts by each other 
# observers
ggplot(data = phys_traits_obs, aes(x = log(Glands), y = log(JM_glands))) +
  geom_abline(slope = 0.98, size = 1, linetype = "dashed") +
  geom_point(aes(color = Obs)) +
  theme_bw()

# Fit alternative models for gland counts.
glands1 <- lm(log(JM_glands) ~ Obs*log(Glands) - 1, data = phys_traits_obs)
glands2 <- lm(log(JM_glands) ~ Obs + log(Glands) - 1, data = phys_traits_obs)
glands3 <- lm(log(JM_glands) ~ Obs:log(Glands) - 1, data = phys_traits_obs)
BIC(glands1, glands2, glands3)

# View best-fit model for gland counts
summary(glands2)
Anova(glands2, type = 3)

# Use coefficient estimates from best-fit model to plot lines by observer
ggplot(data = phys_traits_obs, aes(x = log(Glands), y = log(JM_glands))) +
  geom_abline(slope = 0.98, size = 1, linetype = "dashed") +
  geom_point(aes(color = Obs)) +
  geom_line(aes(x = log(Glands), y = 0.588 + 0.91*log(Glands)),
            color = "red") + #observer JO
  geom_line(aes(x = log(Glands), y = 0.439 + 0.91*log(Glands)),
            color = "forestgreen") + # observer MH
  geom_line(aes(x = log(Glands), y = 0.730 + 0.91*log(Glands)),
            color = "blue") + #observer SB
  theme_bw()

hist(log(phys_traits_obs$Trichomes))

# View relationship between trichome counts by JM and trichome counts by each
# other observer
ggplot(data = phys_traits_obs, aes(x = log(Trichomes), y = log(JM_trichomes))) +
  geom_abline(slope = 1, size = 1, linetype = "dashed") +
  geom_point(aes(color = Obs)) +
  theme_bw()

# Fit alternative models for trichome counts.
trich1 <- lm(log(JM_trichomes) ~ Obs*log(Trichomes) - 1, data = phys_traits_obs)
trich2 <- lm(log(JM_trichomes) ~ Obs + log(Trichomes) - 1, data = phys_traits_obs)
trich3 <- lm(log(JM_trichomes) ~ Obs:log(Trichomes) - 1, data = phys_traits_obs)
BIC(trich1, trich2, trich3)

# View best-fit model for trichome counts.
summary(trich1)
Anova(trich1, type = 3)

# Use coefficient estimates from best-fit model to plot lines by observer
ggplot(data = phys_traits_obs, aes(x = log(Trichomes), y = log(JM_trichomes))) +
  geom_abline(slope = 1, size = 1, linetype = "dashed") +
  geom_point(aes(color = Obs)) +
  geom_line(aes(x = log(Trichomes), y = -0.178 + 1.023*log(Trichomes)),
            color = "red") + #observer JO
  geom_line(aes(x = log(Trichomes), y = 1.351 + (1.023-0.292)*log(Trichomes)),
            color = "forestgreen") + #observer MH
  geom_line(aes(x = log(Trichomes), y = 0.553 + (1.023-0.131)*log(Trichomes)),
            color = "blue") + #observer SB
  theme_bw()


# Read in full leaf physical trait data (trichomes and terpenoid glands) dataset
phys_traits <- read.csv("raw-data/phys-trait-data.csv")
# Remove duplicate measurements and 
phys_traits <- phys_traits[-which(phys_traits$Notes == "duplicate"),]
phys_traits$Glands <- as.numeric(phys_traits$Glands)

# Use coefficient estimates from best-fit models to standardize trichome and 
# gland counts
for (i in 1:length(phys_traits[,1])) {
  if(phys_traits[i, 'Obs'] == "JM") { 
    phys_traits[i, 'Glands_std'] <- phys_traits[i, 'Glands']
    phys_traits[i, 'Trichomes_std'] <- phys_traits[i, 'Trichomes']
    }
  if(phys_traits[i, 'Obs'] == "JO") { 
    phys_traits[i, 'Glands_std'] <- exp(0.588 + 0.91*log(phys_traits[i, 'Glands']))
    phys_traits[i, 'Trichomes_std'] <- exp(-0.178 + 1.023*log(phys_traits[i, 'Trichomes']))
  }
  if(phys_traits[i, 'Obs'] == "MH") {
    phys_traits[i, 'Glands_std'] <- exp(0.439 + 0.91*log(phys_traits[i, 'Glands']))
    phys_traits[i, 'Trichomes_std'] <- exp(1.351 + (1.023-0.292)*log(phys_traits[i, 'Trichomes']))
  }
  if(phys_traits[i, 'Obs'] == "SB") {
    phys_traits[i, 'Glands_std'] <- exp(0.730 + 0.91*log(phys_traits[i, 'Glands']))
    phys_traits[i, 'Trichomes_std'] <- exp(0.553 + (1.023-0.131)*log(phys_traits[i, 'Trichomes']))
  }
}

phys_traits$Glands_std <- as.numeric(phys_traits$Glands_std)

## Convert standardized trichome and gland counts to densities per sq cm, averaged
##  by plant.

# The leaf disc is 5 mm in diameter = 19.63 sq mm. Convert numbers of trichomes
# and terpenoid glands per disc to densities per sq cm by multiplying by 5.09.
phys_traits$Glands_std_cm2 <- with(phys_traits, Glands_std*5.09)
phys_traits$Trichomes_std_cm2 <- with(phys_traits, Trichomes_std*5.09)

# Average trichome and gland densities per sq cm by plant
phys_traits_sum <- doBy::summaryBy(Glands_std_cm2 + Trichomes_std_cm2 ~
                                     Plant_ID, data = phys_traits,
                                   FUN = function(x) { mean(x, na.rm = T) },
                                   keep.names = T)

## Add midvein length data and gland area data to physical traits dataframe.

# Read in midvein lengths for sampled leaves
midvein <- read.csv("raw-data/leaf-trait-midvein-lengths.csv")
midvein_sum <- doBy::summaryBy(Midvein_length ~ Plant_ID,
                               data = midvein, keep.names = T,
                               FUN = function(x) { round(mean(x), digits = 1) })

# Add midvein lengths to physical traits dataframe.
phys_traits_sum$Midvein_length <- midvein_sum[
  match(phys_traits_sum$Plant_ID, midvein_sum$Plant_ID), 'Midvein_length']

# Read in gland area data.
gland_area <- read.csv("raw-data/terpenoid-gland-areas.csv")

# Extract plant ID and image ID from Image_file name.
gland_area$Plant_ID <- with(gland_area, substr(Image_file, start = 12, stop = 15))
gland_area$Image_rep <- with(gland_area, substr(Image_file, start = 19, stop = 19))
gland_area$Plant_ID <- with(gland_area, as.numeric(Plant_ID))

# Average gland areas by plant
gland_area_sum <- doBy::summaryBy(Total_area ~ Plant_ID, data = gland_area,
                                  FUN = function(x) { mean(x, na.rm = T) },
                                  keep.names = T)
# Get prop. of leaf discs covered in glands
gland_area_sum$Gland_prop <- with(gland_area_sum, Total_area/19.63)

# Add gland areas and prop. covered in glands to physical traits dataframe.
phys_traits_sum[,c('Gland_area', 'Gland_prop')] <- gland_area_sum[
  match(phys_traits_sum$Plant_ID, gland_area_sum$Plant_ID),
  c('Total_area', 'Gland_prop')
]

# Save cleaned data as .csv
write.csv(phys_traits_sum, "cleaned-data/phys-trait-std.csv", row.names = FALSE)

rm(phys_traits, phys_traits_obs, phys_traits_sum, gland_area, gland_area_sum, 
   midvein, midvein_sum)
