library(readxl)
library(ggplot2)
library(openxlsx)
library(mvtnorm)
library(survival)
library(MASS)
library(TH.data)
library(tidyr)
library(tidyverse)
library(dplyr)
library(palmerpenguins)
library(emmeans)
library(multcomp)
library(ggsci)
library(tidytext)
library(ggpubr)
library(ggrepel)
library(GGally)
library(PerformanceAnalytics)
library(webr)
library(plotly)
library(stats)
library(plotly)
library(plot3D)
library(plotrix)


# Load required package
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# Set parameters for seasons, stations, and replicates
seasons <- c("Pre-monsoon", "Monsoon", "Post-monsoon")
stations <- c("Station_1", "Station_2")
replicates <- 10

set.seed(123)

# Define function to generate water quality data with specified mean and sd
generate_data <- function(season) {
  if (season == "Pre-monsoon") {
    data.frame(
      Temperature = rnorm(1, mean = runif(1, 29, 33), sd = 1.83),     # Higher temperatures before rains
      Salinity = rnorm(1, mean = runif(1, 6, 9), sd = 1.02),       # Salinity tends to be higher due to evaporation
      pH = rnorm(1, mean = runif(1, 7.0, 8.0), sd = 0.53),         # Slightly alkaline
      DO = rnorm(1, mean = runif(1, 6, 8), sd = 0.32),             # Moderate DO levels
      Transparency = rnorm(1, mean = runif(1, 40, 50), sd = 3.34), # Water is clearer
      PO4P = rnorm(1, mean = runif(1, 0.05, 0.2), sd = 0.08),    # Moderate phosphate levels
      NO2N = rnorm(1, mean = runif(1, 0.02, 0.1), sd = 0.04),    # Low to moderate nitrate
      SiO3Si = rnorm(1, mean = runif(1, 1, 5), sd = 0.5),        # Moderate silicate levels
      TDS = rnorm(1, mean = runif(1, 500, 600), sd = 12.4),         # TDS range based on season
      TSS = rnorm(1, mean = runif(1, 10, 30), sd = 7.34)            # TSS range based on season
    )
  } else if (season == "Monsoon") {
    data.frame(
      Temperature = rnorm(1, mean = runif(1, 25, 30), sd = 2.14),    # Cooler due to rainfall
      Salinity = rnorm(1, mean = runif(1, 2, 6), sd = 0.82),      # Diluted salinity from rainfall
      pH = rnorm(1, mean = runif(1, 6.5, 7.5), sd = 0.53),        # Near neutral due to rainwater input
      DO = rnorm(1, mean = runif(1, 5, 7), sd = 0.55),            # DO may be higher due to runoff and aeration
      Transparency = rnorm(1, mean = runif(1, 30, 40), sd = 5),  # Lower transparency due to runoff
      PO4P = rnorm(1, mean = runif(1, 0.1, 0.3), sd = 0.05),     # Higher phosphates due to runoff
      NO2N = rnorm(1, mean = runif(1, 0.1, 0.2), sd = 0.02),     # Higher nitrates due to runoff
      SiO3Si = rnorm(1, mean = runif(1, 3, 10), sd = 0.5),       # Increased silicates
      TDS = rnorm(1, mean = runif(1, 200, 300), sd = 15.43),         # TDS lowers due to dilution
      TSS = rnorm(1, mean = runif(1, 50, 200), sd = 8.54)           # TSS increases due to runoff
    )
  } else { # Post-monsoon
    data.frame(
      Temperature = rnorm(1, mean = runif(1, 26, 31), sd = 1.6),    # Moderate temperatures
      Salinity = rnorm(1, mean = runif(1, 1, 4), sd = 0.98),      # Salinity stabilizes
      pH = rnorm(1, mean = runif(1, 7.0, 8.2), sd = 0.34),        # Alkaline after sediment settling
      DO = rnorm(1, mean = runif(1, 6, 9), sd = 0.55),            # DO stabilizes
      Transparency = rnorm(1, mean = runif(1, 50, 70), sd = 6.2), # Water becomes clearer again
      PO4P = rnorm(1, mean = runif(1, 0.05, 0.15), sd = 0.04),   # Phosphate decreases as runoff declines
      NO2N = rnorm(1, mean = runif(1, 0.02, 0.1), sd = 0.03),    # Nitrate decreases
      SiO3Si = rnorm(1, mean = runif(1, 1, 5), sd = 0.56),        # Silicate levels moderate
      TDS = rnorm(1, mean = runif(1, 300, 500), sd = 13.65),        # TDS stabilizes
      TSS = rnorm(1, mean = runif(1, 20, 100), sd = 6.43)          # TSS decreases as sediment settles
    )
  }
}


# Generate data for each season, station, and replicate
water_quality_data <- expand.grid(Season = seasons, Station = stations, Replicate = 1:replicates) %>%
  rowwise() %>%
  mutate(data = list(generate_data(Season))) %>%
  unnest(data)

# Write data to CSV file
write.csv(water_quality_data, "Sangu_River_Water_Quality_Data.csv", row.names = FALSE)


# Generate data for each season, station, and replicate
water_quality_data <- expand.grid(Season = seasons, Station = stations, Replicate = 1:replicates) %>%
  rowwise() %>%
  mutate(data = list(generate_data(Season))) %>%
  unnest(data)

# Write data to CSV file
write.csv(water_quality_data, "Sangu_River_Water_Quality_Data.csv", row.names = FALSE)



summary3 <- read_excel("ZooSangu.xlsx", sheet = "water")


sum <-   summary3 %>% 
  pivot_longer(cols = 2:11, names_to = "Parameters", values_to = "value" ) %>% 
  group_by(Season, Parameters) %>% 
  summarize(mean = round(mean(value), digits = 2), sd = round(sd(value), digits = 2), .groups = "drop") %>% 
  mutate(mean_pm = paste0( mean, "±", sd))
 

#
write.xlsx(sum, "summarywq.xlsx")


# let's dive into the analysis
zoo <- read_excel("ZooSangu.xlsx" , sheet = "water")

# all boxplot
zoo %>% 
  pivot_longer(cols = 2:11, names_to = "Parameters", values_to = "values") %>% 
  ggplot(aes(Season, values, fill = Season))+
  geom_boxplot(alpha = 1, show.legend = T, outlier.shape = NA)+
  geom_jitter(size = 1, width = 0.2, alpha = 0.5, color = "black", show.legend = T)+
  geom_violin( alpha = 0.2, color = "black", show.legend = T)+
  facet_wrap(~Parameters, scales = "free_y")+
  stat_summary(fun = mean, geom = "point", size = 3, color = "yellow", show.legend = T)+
  theme_bw()

  #scale_fill_brewer(palette = "Set2") 
  
# for individual boxplot i'm making a long data set for all 

ldata <- zoo %>% 
  pivot_longer(cols = 2:11, names_to = "Parameters", values_to = "values")
  

# for Temp
adata <- ldata %>% 
  filter(Parameters == "Temp") %>% 
  select(1,3)

OneWayAnova(data = adata)


# for Salinity
adata <- ldata %>% 
  filter(Parameters == "Salinity") %>% 
  select(1,3)

OneWayAnova(data = adata)


# for pH           
adata <- ldata %>% 
  filter(Parameters == "pH") %>% 
  select(1,3)

OneWayAnova(data = adata)


# for DO           
adata <- ldata %>% 
  filter(Parameters == "DO") %>% 
  select(1,3)

OneWayAnova(data = adata)

# for Trans           
adata <- ldata %>% 
  filter(Parameters == "Trans") %>% 
  select(1,3)

OneWayAnova(data = adata)

# for PO4P           
adata <- ldata %>% 
  filter(Parameters == "PO4P") %>% 
  select(1,3)

OneWayAnova(data = adata)

# for NO2N           
adata <- ldata %>% 
  filter(Parameters == "NO2N") %>% 
  select(1,3)

OneWayAnova(data = adata)


# for SiO3Si           
adata <- ldata %>% 
  filter(Parameters == "SiO3Si") %>% 
  select(1,3)

OneWayAnova(data = adata)

# for TDS           
adata <- ldata %>% 
  filter(Parameters == "TDS") %>% 
  select(1,3)

OneWayAnova(data = adata)

# for TSS           
adata <- ldata %>% 
  filter(Parameters == "TSS") %>% 
  select(1,3)

OneWayAnova(data = adata)


# abundance data

abundance <- read_excel("ZooSangu.xlsx", sheet = 1)

azdata <- abundance %>% 
  pivot_longer(cols = 2:4, names_to = "Seasons" , values_to = "IndPm3" )
  

azdata %>% 
  ggplot(aes(Seasons, IndPm3, fill = Species))+
  geom_bar(stat = "identity", position = "stack")+
  theme_bw()+
  scale_fill_npg()

azdata %>% 
  ggplot(aes(Seasons, IndPm3, fill = Species))+
  geom_bar(stat = "identity", position = "fill")+
  theme_bw()+
  scale_fill_npg()




#PCA of ZooSangu

zpca <- read_excel("ZooSangu.xlsx", sheet = "PCA")

#mydata <- read_excel("Exam.xlsx", sheet = "PCA")
pcadata <- zpca %>% 
  select(2:ncol(zpca))


#Data standardization
#Note that, by default, the function PCA() [in FactoMineR], standardizes the data automatically 
#during the PCA; so you don’t need do this transformation before the PCA.

#The function PCA() [FactoMineR package] can be used. A simplified format is :
#PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE)
library("FactoMineR")
res.pca <- PCA(pcadata, ncp = 5, graph = F)

#Eigenvalues / Variances
library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#Graph of variables

var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)


#Correlation circle

# Coordinates of variables
head(var$coord, 4)

fviz_pca_var(res.pca, col.var = "black", axes = c(1 , 2), repel = T)


#Quality of representation
head(var$cos2, 4)
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

#variables with low cos2 values will be colored in “white”
#variables with mid cos2 values will be colored in “blue”
#variables with high cos2 values will be colored in red


# Color by cos2 values: quality on the factor map for axis 1 and 2
fviz_pca_var(res.pca, axes = c(1,2), col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)


# Color by cos2 values: quality on the factor map for axis 1 and 3
fviz_pca_var(res.pca, axes = c(1,3), col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Change the transparency by cos2 values
fviz_pca_var(res.pca, alpha.var = "cos2")


#



# for diversity


# diversity

zdata <- read_excel("ZooSangu.xlsx", sheet = "Diversity")

# all boxplot
zdata %>% 
  #pivot_longer(cols = 2:10, names_to = "Seasons", values_to = "values") %>% 
  ggplot(aes(Seasons, values, fill =  Seasons))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, width = 0.2, alpha = 0.5, color = "black", show.legend = F)+
  geom_violin( alpha = 0.2, color = "black", show.legend = F)+
  facet_wrap(~Index, scales = "free_y")+
  stat_summary(fun = mean, geom = "point", size = 3, color = "yellow", show.legend = F)+
  theme_bw()




# for Simpson   
adata <- zdata %>% 
  filter(Index  == "Simpson") %>% 
  select(1,3)

OneWayAnova(data = adata)



# for Simpson   
adata <- zdata %>% 
  filter(Index  == "Simpson") %>% 
  select(1,3)

OneWayAnova(data = adata)

Shannon


# for Margalef   
adata <- zdata %>% 
  filter(Index  == "Margalef") %>% 
  select(1,3)
OneWayAnova(data = adata)

# for Evenness   
adata <- zdata %>% 
  filter(Index  == "Evenness") %>% 
  select(1,3)
OneWayAnova(data = adata)


# for Shannon   
adata <- zdata %>% 
  filter(Index  == "Shannon") %>% 
  select(1,3)
OneWayAnova(data = adata)

# nmds

library(vegan)
ndata <- read_excel("ZooSangu.xlsx", sheet = "nmds")

perm_data <- ndata[, 2:ncol(ndata)]

# Square root Transformation


perm_data_trans <- sqrt(perm_data)


write.xlsx(perm_data_trans, "sqrt.xlsx" )
#Rub NMDS model for visualization the composition

nmds_result <- metaMDS(perm_data_trans, k = 2, distance = "bray", trymax = 10000 )

nmds_result$stress

plot(nmds_result)


#Extracting the NMDS scores
nmds_scores <- as.data.frame(scores(nmds_result)$sites)

ggplot(nmds_scores, aes(NMDS1, NMDS2,shape = ndata$Seasons, col = row.names(perm_data)))+
  geom_point(size = 4, show.legend = F)+
  geom_text(label = ndata$Seasons, nudge_x = 0.015, nudge_y = 0.015, show.legend = F)+
  #theme(legend.position = "")+
  theme_bw()


# db-rda analysis

library(vegan)
myrdata <- read_excel("ZooSangu.xlsx", sheet = "RDA")

tdata <- myrdata %>% 
  select(2:ncol(myrdata)) %>% 
  mutate(across( -pH, ~log10(.+1))) #%>% 
View()

spdata <- tdata %>% 
  select(1:11) %>% 
  as.matrix()

rownames(spdata) <- c("Pre",
                      "mon",
                      "Post")
pdata <-  as.data.frame(spdata)


endata <- tdata %>% 
  select(12:ncol(tdata)) %>% 
  as.matrix()


rownames(endata) <- c("Pre",
                      "mon",
                      "Post")
edata <- as.data.frame(endata)

seasons <- factor(c("Premonsoon", "Monsoon", "PostMonsoon"))
rda_result <- rda(pdata ~ ., data = edata)

plot(rda_result, scaling = 2) # or scaling = 1

vif.cca(rda_result)

# Summary of db-RDA result
summary(rda_result)
# Extract scores for species and environmental variables
species_scores <- as.data.frame(scores(rda_result, display = "species", scaling = 2 ))
env_scores <- as.data.frame(scores(rda_result, display = "bp", scaling = 2))
sample_scores <- as.data.frame(scores(rda_result, display = "sites", scaling = 2))

# Add season information to sample scores
sample_scores$Season <- seasons

# Extract eigenvalues for axis percentages
eigenvalues <- summary(rda_result)$cont$importance[2,]
axis1_percent <- round(eigenvalues[1] * 100, 2)
axis2_percent <- round(eigenvalues[2] * 100, 2)
# Plot using ggplot2
ggplot() +
  geom_point(data = sample_scores, aes(x = RDA1, y = RDA2, shape = Season, color = Season), size = 3) +
  geom_segment(data = env_scores, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "blue") +
  geom_segment(data = species_scores, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text_repel(data = env_scores, aes(x = RDA1, y = RDA2, label = rownames(env_scores)), 
                  color = "blue", vjust = -1) +
  geom_text_repel(data = species_scores, aes(x = RDA1, y = RDA2, label = rownames(species_scores)), 
                  color = "red", vjust = 1) +
  labs( 
    x = paste0("RDA 1 (", axis1_percent, "%)"), 
    y = paste0("RDA 2 (", axis2_percent, "%)")) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")




# correlation plot in r
