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


#winter pH
data <- penguins %>% 
  select(species, bill_length_mm)



OneWayAnova <- function(data){
  aov <- aov(data[[2]] ~ data[[1]]) # doing anova with the data frame
  summary <- summary(aov) # summary of the anova model
  TukeySHD <- TukeyHSD(aov) # tukey shd test
  tukyplot <- plot(TukeySHD) # ploting the tukey result
  emmean <- emmeans::emmeans(aov, specs = "data[[1]]")
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  
  plot <- ggplot2::ggplot(data = data, aes(x = as.factor(data[[1]]), y = data[[2]], fill = as.factor(data[[1]]))) +
    geom_boxplot(show.legend = F, outlier.shape = NA, alpha = 1) +
   # geom_jitter(size = 3, show.legend = F, alpha = 0.6) +
    theme_bw() +
    # Increase text size for geom_text
    geom_text(data = emmean_cld, aes(x = as.factor(emmean_cld[[1]]), y = upper.CL, label = .group), 
              vjust = -0.5, size = 8, inherit.aes = F) # Adjust 'size' for geom_text
  
  return(list(summary, emmean_cld, plot))
    #list(summary, TukeySHD, emmean, emmean_cld, plot)
}
myzoo <- read_excel("Zoo1.xlsx", sheet = "water")

mydata <- myzoo %>% 
  select(2,3)

OneWayAnova(data = myzoo)

data <- wdata %>% 
  filter(Parameters == "pH") %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "values") %>% 
  select(Seasons, Stations, values)



TwoWayAnova <- function(data){
  aov <- aov(data[[3]] ~ data[[2]] + data[[1]]) # doing anova with the data frame
  summary <- summary(aov) # summary of the anova model
  TukeySHD <- TukeyHSD(aov) # tukey shd test
  tukyplot <- plot(TukeySHD) # ploting the tukey result
  emmean <- emmeans::emmeans(aov, specs = c("data[[2]]", "data[[1]]"))
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  
  plot <-  ggplot(data = data, aes(x = as.factor(data[[2]]), y = data[[3]], fill = as.factor(data[[1]]))) +
    geom_boxplot(show.legend = F, outlier.shape = NA, alpha = 1) +
    # geom_jitter(size = 3, show.legend = F, alpha = 0.6) +
    theme_bw() +
    # Increase text size for geom_text
    geom_text(data = emmean_cld, aes(x = as.factor(emmean_cld[[1]]), y = upper.CL, label = .group),
              position = position_dodge(width = 0.75) ,vjust = -0.5, size = 8, inherit.aes = F) # Adjust 'size' for geom_text
  

  return(list(summary, emmean_cld, plot))
  #list(summary, TukeySHD, emmean, emmean_cld, plot)
}


TwoWayAnova(data = dataz)

library(aLBI)

lfq <- lenfreq01

FishPar(data = lfq, resample = 1000, progress = F)

# start

zdata <- read_excel("Zoo1.xlsx", sheet = "water")

wsum <- zdata %>% 
  pivot_longer(cols = 3:4, 
               names_to = "Seasons",
               values_to = "Values") %>% 
  group_by(Stations, Parameters, Seasons) %>% 
  summarise(mean = mean(Values), sd = sd(Values))

write.xlsx(wsum, "wsum.xlsx")

ws <- zdata %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "values") %>% 
  group_by(Stations, Seasons, Parameters) %>% 
  summarise(mean = mean(values)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Parameters, values_from = mean)

write.xlsx(ws,"ws.xlsx",  overwrite = T)


zdata %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "Value") %>% 
  ggplot(aes(Stations, Value, fill = Seasons))+
  geom_boxplot(show.legend = T)+
  facet_wrap(~Parameters, scales = "free_y")+
  theme_bw()+
  theme(legend.position = "top")




# for pH
adata <- zdata %>% 
  filter(Parameters == "pH") %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "values")

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)

# for temp
adata <- zdata %>% 
  filter(Parameters == "Temp") %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "values")

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


# for salinity
adata <- zdata %>% 
  filter(Parameters == "Salinity") %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "values")

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)



# for DO
adata <- zdata %>% 
  filter(Parameters == "DO") %>% 
  pivot_longer(cols = 3:4, names_to = "Seasons", values_to = "values")

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)



# abundance data

azdata <- read_excel("Zoo1.xlsx", sheet = 3)

z <- azdata %>% 
  select(1:3, 6) %>% 
  pivot_wider(names_from = Species, values_from = rIndPm)
           
wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, "sheet1", z)
addWorksheet(wb, "sheet2")
writeData(wb, "sheet2", ws)

saveWorkbook(wb, "zoow.xlsx", overwrite = T)

azdata %>% 
  ggplot(aes(Stations, IndPm3, fill = Species))+
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(~Seasons)+
  theme_bw()+
  scale_fill_npg()

azdata %>% 
  ggplot(aes(Stations, IndPm3, fill = Species))+
  geom_bar(stat = "identity", position = "fill")+
  facet_wrap(~Seasons)+
  theme_bw()+
  scale_fill_npg()




#PCA of zoow

zpca <- read_excel("zoow.xlsx", sheet = 1)

#mydata <- read_excel("Exam.xlsx", sheet = "PCA")
pcadata <- zpca %>% 
  select(3:ncol(zpca))


#Data standardization
#Note that, by default, the function PCA() [in FactoMineR], standardizes the data automatically 
#during the PCA; so you don’t need do this transformation before the PCA.

#The function PCA() [FactoMineR package] can be used. A simplified format is :
#PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE)
library("FactoMineR")
res.pca <- PCA(pcadata, graph = F)

#Visualization and Interpretation

#get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
#fviz_eig(res.pca): Visualize the eigenvalues
#get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
#fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
#fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.



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

fviz_pca_var(res.pca, col.var = "black")

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


# corelation

library(PerformanceAnalytics)

library("PerformanceAnalytics")
coor <- blg %>% 
  select(2:5)

chart.Correlation(pcadata, histogram = F, pch= 19)


devtools::install_github("Hy4m/linkET", force = T)

library(linkET)


# for diversity
dv <- zpca %>% 
  select(1:13) %>% 
  pivot_longer(cols = 3:13, names_to = "Species", values_to = "Number") %>% 
  pivot_wider(names_from = Stations, values_from = Number)
  
  
write.xlsx(dv, "dv.xlsx", overwrite = T)  



# diversity

zdata <- read_excel("dv.xlsx", sheet = 3)

sdiv <- zdata %>% 
  group_by(Seasons, Index, Stations) %>% 
  summarize(mean = mean(values),
            sd = sd(values))

write.xlsx(sdiv, "divs.xlsx")

zdata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = T)+
  facet_wrap(~Index, scales = "free_y")+
  theme_bw()+
  theme(legend.position = "top")




# for Simpson   
adata <- zdata %>% 
  filter(Index  == "Simpson")

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)

# for temp
adata <- zdata %>% 
  filter(Index == "Shannon_H") 

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


# for salinity
adata <- zdata %>% 
  filter(Index == "Evenness") 

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)



# for Margalef  
adata <- zdata %>% 
  filter(Index == "Margalef") 

#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
emmean_cld

# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)





# nmds

ndata <- read_excel("zoow.xlsx", sheet = 1)

perm_data <- ndata[, 3:14]

# Square root Transformation


perm_data_trans <- sqrt(perm_data)

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


# rda analysis

library(vegan)

myrdata <- read_excel("zoow.xlsx", sheet = 2)

tdata <- myrdata %>% 
  select(2:ncol(myrdata)) %>% 
  mutate(across( -pH, ~log10(.+1))) #%>% 
View()

spdata <- tdata %>% 
  select(1:11) %>% 
  as.matrix()

rownames(spdata) <- c("Pre1", "Pre2", "Pre3",
                      "Post1", "Post2", "Post3")
pdata <-  as.data.frame(spdata)


endata <- tdata %>% 
  select(12:ncol(tdata)) %>% 
  as.matrix()


rownames(endata) <- c("Pre1", "Pre2", "Pre3",
                      "Post1", "Post2", "Post3" )
edata <- as.data.frame(endata)

seasons <- factor(c("Premonsoon", "PostMonsoon"))
rda_result <- rda(pdata ~ ., data = edata)

plot(rda_result)
# Summary of RDA result
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

library(ggcorrplot)

ndata

corr <- ndata %>% 
  select(3:ncol(ndata)) %>% as.matrix()
sdata <- ndata %>% 
  select(where(is.numeric))

# Load necessary libraries
library(corrplot)
library(dplyr)

# Step 1: Select only numeric columns (remove 'Stations' and 'Seasons')
numeric_sdata <- sdata %>% select(where(is.numeric))

# Step 2: Compute the correlation matrix
cor_matrix <- cor(numeric_sdata, use = "complete.obs")

# Step 3: Visualize the correlation matrix using corrplot
corrplot(cor_matrix,
         method = "circle", 
         type = "lower",
         tl.col = "black", 
         tl.srt = 45)

