library(vegan)
library(readxl)



# Principal component analysis --------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyverse)
library(usethis)
library(devtools)

#install.packages("ggbiplot")
head(iris)

#install_github("vqv/ggbiplot")

pc <-  iris %>%
  select(1:4) %>% 
  prcomp(center = T, scale = T)

pc$scale 

summary(pc)

ggbiplot(pc, 
         obs.scale = 1,
         var.scale = 1, 
         groups = iris$Species,
         ellipse = T,
         circle = T,
         ellipse.prob = 0.68
)+
  theme(panel.background = element_blank(),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black")
        
  )



#nMDS analysis

library(vegan)
library(ggplot2)
library(tidyverse)


set.seed(100)

nmds <- read.csv("Data.csv")

nmds[, 2:7]
ndata <- nmds %>% 
  select(2:7)

nmds_result <- metaMDS(ndata, distance = "bray", k = 2, trymax = 100)
plot(nmds_result)

#Stress plot

stressplot(nmds_result)
nmds_result$stress  
# stress value > 0.2 indicates a poor fit
# stress value 0.05 - 0.1 indicates a good fit
# stress value < 0.05 indicates excellent fit

# ploting using ggplot2

# with my data

library(vegan)
library(ggplot2)
library(dendextend)
library(tidyverse)
library(dplyr)

# Sample data: abundance of genera in different samples
# Load required libraries



#PREMANOVA in r
#Load paclages

library(vegan)
library(tidyverse)
library(readxl)

#import dataset

set.seed(100)

data <- read.csv("Dataset.csv")

data <- read_excel("Phyto.xlsx", sheet = "p")

perm_data <- data[, 3:14]

# Square root Transformation


perm_data_trans <- sqrt(perm_data)

#Rub NMDS model for visualization the composition

nmds_result <- metaMDS(perm_data_trans, k = 2, distance = "bray", trymax = 10000 )

nmds_result$stress

plot(nmds_result)


#Extracting the NMDS scores
nmds_scores <- as.data.frame(scores(nmds_result)$sites)

ggplot(nmds_scores, aes(NMDS1, NMDS2,shape = data$Seasons, col = row.names(perm_data)))+
  geom_point(size = 4, show.legend = F)+
  geom_text(label = data$Seasons, nudge_x = 0.015, nudge_y = 0.015, show.legend = F)+
  #theme(legend.position = "")+
  theme_bw()



#PCA 
#install.packages("aplot")
  
  
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
theme_set(theme_bw())
data(dataset)
# PCoA

t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
# extract the axis scores
tmp <- t1$res_ordination$scores
# differential test with trans_env class
t2 <- trans_env$new(dataset = dataset, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "Group", method = "anova")
  
p1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
p1
# groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
p2 <- t2$plot_diff(measure = "PCo1", add_sig = T) + theme_bw() + coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2
 p3 <- t2$plot_diff(measure = "PCo2", add_sig = T) + theme_bw() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
p3


#Differentiate test
library(microeco)
data(dataset)
t1 <- trans_func$new(dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# use list to prepare data
tmp <- list()
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp$func <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the list as taxa_abund in your microtable object
dataset$taxa_abund <- tmp
# use trans_diff class to perform differential test
t2 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "all")
t2$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)")



library(microeco)
data(dataset)
data(env_data_16S)
t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
t1$cal_ordination(method = "RDA", taxa_level = "Genus")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
t1$trans_ordination(adjust_arrow_length = TRUE)
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
g1
#ggplot2::ggsave("RDA.pdf", g1, width = 8, height = 6.5)
# use capture.output to save output
capture.output(t1$res_ordination_R2, file = "RDA_R2.txt")
capture.output(t1$res_ordination_envfit, file = "RDA_envfit.txt")
# save data.frame objects
write.table(t1$res_ordination_terms, "RDA_anova_termsig.txt", sep = "\t")
write.table(t1$res_ordination_axis, "RDA_anova_axissig.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_sites, "RDA_axis_sample.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_arrows, "RDA_axis_term.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_arrows_spe, "RDA_axis_taxa.txt", sep = "\t")






#PCA 
#install.packages("FactoMineR")
#install.packages("factoextra")


library("FactoMineR")
library("factoextra")
data(decathlon2)
head(decathlon2)
dim(decathlon2)
View(decathlon2)

mydata <- read_excel("Exam.xlsx", sheet = "PCA")
pcadata <- mydata %>% 
  select(2:ncol(mydata))
 dim(pcadata) 

decathlon2.active <- decathlon2[1:23, 1:10]
dim(decathlon2.active)
View(decathlon2.active)
head(decathlon2.active[, 1:6], 4)

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


#Contributions of variables to PCs

head(var$contrib, 4)
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)  

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#The total contribution to PC1 and PC2 is obtained with the following R code:
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)

#The most important (or, contributing) variables can be highlighted on the correlation plot as follow:

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T
)


## Change the transparency by contrib values
fviz_pca_var(res.pca, alpha.var = "contrib")


#Color by a custom continuous variable

# Create a random continuous variable of length 10
set.seed(123)
my.cont.var <- rnorm(16)
# Color variables by the continuous variable
fviz_pca_var(res.pca, col.var = my.cont.var,
             gradient.cols = c("blue", "yellow", "red"),
             legend.title = "Cont.Var",
            repel = T  )

#Color by groups

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")


#Dimension description
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1

# Description of dimension 2
res.desc$Dim.2
#Graph of individuals
ind <- get_pca_ind(res.pca)
ind

#to get access to the different components use this:
# Coordinates of individuals
head(ind$coord)
# Quality of individuals
head(ind$cos2)
# Contributions of individuals
head(ind$contrib)

#Plots: quality and contribution
#The fviz_pca_ind() is used to produce the graph of individuals. To create a simple plot, type this:
  
fviz_pca_ind(res.pca)

#Like variables, it’s also possible to color individuals by their cos2 values:
  
  fviz_pca_ind(res.pca, col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE # Avoid text overlapping (slow if many points)
  )

#Note that, individuals that are similar are grouped together on the plot.

#You can also change the point size according the cos2 of the corresponding individuals:
  
fviz_pca_ind(res.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#To change both point size and color by cos2, try this:
  
fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#To visualize the contribution of individuals to the first two principal components, type this:
  
# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

#Color by a custom continuous variable
# Create a random continuous variable of length 23,
# Same length as the number of active individuals in the PCA
set.seed(123)
my.cont.var <- rnorm(16)
# Color individuals by the continuous variable
fviz_pca_ind(res.pca, col.ind = my.cont.var,
             gradient.cols = c("blue", "yellow", "red"),
             legend.title = "Cont.Var")


#Color by groups
head(iris, 3)
# The variable Species (index = 5) is removed
# before PCA analysis
iris.pca <- PCA(iris[,-5], graph = FALSE)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


# with my data as group
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = mydata$Seasons, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
#To remove the group mean point, specify the argument mean.point = FALSE.

#If you want confidence ellipses instead of concentration ellipses, use ellipse.type = “confidence”.

# Add confidence ellipses
fviz_pca_ind(iris.pca, geom.ind = "point", col.ind = iris$Species, 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Groups"
)


# with my data
fviz_pca_ind(res.pca, geom.ind = "point", col.ind = mydata$Seasons, 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Groups"
)
#Note that, allowed values for palette include:
  
#“grey” for grey color palettes;
#brewer palettes e.g. “RdBu”, “Blues”, …; To view all, type this in R: RColorBrewer::display.brewer.all().
#custom color palette e.g. c(“blue”, “red”);
#and scientific journal palettes from ggsci R package, e.g.: “npg”, “aaas”, “lancet”, “jco”, “ucscgb”, “uchicago”, “simpsons” and “rickandmorty”.
#For example, to use the jco (journal of clinical oncology) color palette, type this:

fviz_pca_ind(iris.pca,
             label = "none", # hide individual labels
             habillage = iris$Species, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = "jco"
)

#Graph customization

#Note that, fviz_pca_ind() and fviz_pca_var() and related functions are 
#wrapper around the core function fviz() [in factoextra]. fviz() is a wrapper 
#around the function ggscatter() [in ggpubr]. Therefore, further arguments, 
#to be passed to the function fviz() and ggscatter(), can be specified in 
#fviz_pca_ind() and fviz_pca_var().

#Dimension
# Variables on dimensions 2 and 3
fviz_pca_var(res.pca, axes = c(2, 3))
# Individuals on dimensions 2 and 3
fviz_pca_ind(res.pca, axes = c(2, 3))

#Plot elements: point, text, arrow
#geom.var: a text specifying the geometry to be used for plotting variables. Allowed values are the combination of c(“point”, “arrow”, “text”).
#Use geom.var = "point", to show only points;
#Use geom.var = "text" to show only text labels;
#Use geom.var = c("point", "text") to show both points and text labels
#Use geom.var = c("arrow", "text") to show arrows and labels (default).

# Show variable points and text labels
fviz_pca_var(res.pca, geom.var = c("point", "text"))

#geom.ind: a text specifying the geometry to be used for plotting individuals. Allowed values are the combination of c(“point”, “text”).
#Use geom.ind = "point", to show only points;
#Use geom.ind = "text" to show only text labels;
#Use geom.ind = c("point", "text") to show both point and text labels (default)
# Show individuals text labels only
fviz_pca_ind(res.pca, geom.ind =  "text")

#Size and shape of plot elements
#labelsize: font size for the text labels, e.g.: labelsize = 4.
#pointsize: the size of points, e.g.: pointsize = 1.5.
#arrowsize: the size of arrows. Controls the thickness of arrows, e.g.: arrowsize = 0.5.
#pointshape: the shape of points, pointshape = 21. Type ggpubr::show_point_shapes() to see available point shapes.
# Change the size of arrows an labels
fviz_pca_var(res.pca, arrowsize = 1, labelsize = 5, 
             repel = TRUE)
# Change points size, shape and fill color
# Change labelsize
fviz_pca_ind(res.pca, 
             pointsize = 3, pointshape = 21, fill = "lightblue",
             labelsize = 5, repel = TRUE)

#Ellipses
#As we described in the previous section @ref(color-ind-by-groups), 
#when coloring individuals by groups, you can add point concentration ellipses 
#using the argument addEllipses = TRUE.

#Note that, the argument ellipse.type can be used to change the type of ellipses. Possible values are:
  
#"convex": plot convex hull of a set o points.
#"confidence": plot confidence ellipses around group mean points as the function coord.ellipse() [in FactoMineR].
#"t": assumes a multivariate t-distribution.
#"norm": assumes a multivariate normal distribution.
#"euclid": draws a circle with the radius equal to level, representing the euclidean distance from the center. This ellipse probably won’t appear circular unless coord_fixed() is applied.
#The argument ellipse.level is also available to change the size of the concentration ellipse in normal probability. For example, specify ellipse.level = 0.95 or ellipse.level = 0.66.

# Add confidence ellipses
fviz_pca_ind(iris.pca, geom.ind = "point", 
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Groups"
)
# Convex hull
fviz_pca_ind(iris.pca, geom.ind = "point",
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "convex",
             legend.title = "Groups"
)

#Group mean points
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (but not "text")
             group.ind = iris$Species, # color by groups
             legend.title = "Groups",
             mean.point = FALSE)

#Axis lines
#The argument axes.linetype can be used to specify the line type of axes.
#Default is “dashed”. Allowed values include “blank”, “solid”, “dotted”, etc. 
#To see all possible values type ggpubr::show_line_types() in R.

#To remove axis lines, use axes.linetype = “blank”:
fviz_pca_var(res.pca, axes.linetype = "blank")


#Graphical parameters
#To change easily the graphical of any ggplots, you can use the function ggpar() [ggpubr package]

#The graphical parameters that can be changed using ggpar() include:
  
#Main titles, axis labels and legend titles
#Legend position. Possible values: “top”, “bottom”, “left”, “right”, “none”.
#Color palette.
#Themes. Allowed values include: theme_gray(), theme_bw(), theme_minimal(), theme_classic(), theme_void().

ind.p <- fviz_pca_ind(iris.pca, geom = "point", col.ind = iris$Species)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Iris data set",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Species", legend.position = "top",
              ggtheme = theme_gray(), palette = "jco"
)

#Biplot
#To make a simple biplot of individuals and variables, type this:
  
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969" ,
                repele = T # Individuals color
)

#Note that, the biplot might be only useful when there is a low number of variables and individuals in the data set; otherwise the final plot would be unreadable.

#Note also that, the coordinate of individuals and variables are not constructed on the same space. Therefore, in the biplot, you should mainly focus on the direction of variables but not on their absolute positions on the plot.

#Roughly speaking a biplot can be interpreted as follow:
  
#an individual that is on the same side of a given variable has a high value for this variable;
#an individual that is on the opposite side of a given variable has a low value for this variable.
#Now, using the iris.pca output, let’s :
  
#make a biplot of individuals and variables
#change the color of individuals by groups: col.ind = iris$Species
#show only the labels for variables: label = "var" or use geom.ind = "point"
fviz_pca_biplot(iris.pca, 
                col.ind = iris$Species, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species") 

#To customize individuals and variable colors, we use the helper 
#functions fill_palette() and color_palette() [in ggpubr package].

fviz_pca_biplot(iris.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = iris$Species,
                col.ind = "black",
                # Color variable by groups
                col.var = factor(c("sepal", "sepal", "petal", "petal")),
                
                legend.title = list(fill = "Species", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


#Another complex example is to color individuals by groups (discrete color)
#and variables by their contributions to the principal components (gradient colors). 
#Additionally, we’ll change the transparency of variables by their 
#contributions using the argument alpha.var.
fviz_pca_biplot(iris.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = iris$Species, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Species", color = "Contrib",
                                    alpha = "Contrib")
)


#Supplementary elements
#Specification in PCA
PCA(X, ind.sup = NULL, 
    quanti.sup = NULL, quali.sup = NULL, graph = TRUE)

#X : a data frame. Rows are individuals and columns are numeric variables.
#ind.sup : a numeric vector specifying the indexes of the supplementary individuals
#quanti.sup, quali.sup : a numeric vector specifying, respectively, the indexes of the quantitative and qualitative variables
#graph : a logical value. If TRUE a graph is displayed.

res.pca <- PCA(decathlon2, ind.sup = 24:27, 
               quanti.sup = 11:12, quali.sup = 13, graph=FALSE)


#Quantitative variables
#Predicted results (coordinates, correlation and cos2) for the supplementary quantitative variables:
res.pca$quanti.sup

#Visualize all variables (active and supplementary ones):
fviz_pca_var(res.pca)

#Note that, by default, supplementary quantitative variables are shown in blue color and dashed lines.
#Further arguments to customize the plot:
  
# Change color of variables
fviz_pca_var(res.pca,
             col.var = "black",     # Active variables
             col.quanti.sup = "red" # Suppl. quantitative variables
)
# Hide active variables on the plot, 
# show only supplementary variables
fviz_pca_var(res.pca, invisible = "var")
# Hide supplementary variables
fviz_pca_var(res.pca, invisible = "quanti.sup")

# Plot of active variables
p <- fviz_pca_var(res.pca, invisible = "quanti.sup")
# Add supplementary active variables
fviz_add(p, res.pca$quanti.sup$coord, 
         geom = c("arrow", "text"), 
         color = "red")

#Individuals
#Predicted results for the supplementary individuals (ind.sup):
res.pca$ind.sup
p <- fviz_pca_ind(res.pca, col.ind.sup = "blue", repel = TRUE)
p <- fviz_add(p, res.pca$quali.sup$coord, color = "red")
p


#Qualitative variables
res.pca$quali

fviz_pca_ind(res.pca, habillage = 13,
             addEllipses =TRUE, ellipse.type = "confidence",
             palette = "jco", repel = TRUE) 

#Filtering results
#If you have many individuals/variable, it’s possible to visualize only
#some of them using the arguments select.ind and select.var.

#select.ind, select.var: a selection of individuals/variable to be 
#plotted. Allowed values are NULL or a list containing the arguments name, cos2 or contrib:
  
#name: is a character vector containing individuals/variable names to be plotted
#cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are plotted
#if cos2 > 1, ex: 5, then the top 5 active individuals/variables and top 5 supplementary columns/rows with the highest cos2 are plotted
#contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contributions are plotted
# Visualize variable with cos2 >= 0.6
fviz_pca_var(res.pca, select.var = list(cos2 = 0.6))
# Top 5 active variables with the highest cos2
fviz_pca_var(res.pca, select.var= list(cos2 = 5))
# Select by names
name <- list(name = c("Long.jump", "High.jump", "X100m"))
fviz_pca_var(res.pca, select.var = name)
# top 5 contributing individuals and variable
fviz_pca_biplot(res.pca, select.ind = list(contrib = 5), 
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())


# Print the plot to a pdf file
pdf("myplot.pdf")
print(myplot)
dev.off()

# Scree plot
scree.plot <- fviz_eig(res.pca)
# Plot of individuals
ind.plot <- fviz_pca_ind(res.pca)
# Plot of variables
var.plot <- fviz_pca_var(res.pca)

pdf("PCA.pdf") # Create a new pdf device
print(scree.plot)
print(ind.plot)
print(var.plot)
dev.off() # Close the pdf device

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/


#CCA

#Graph of contingency tables and chi-square test
#install.packages("gplots")
library("gplots")
# 1. convert the data as a table
dt <- as.table(as.matrix(housetasks))
# 2. Graph
balloonplot(t(dt), main ="housetasks", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)
chisq <- chisq.test(dt)
summary(chisq)

mycadata <- read_excel("Phyto.xlsx", sheet = "CA")

mca <- mydata %>% 
  select(3:18) %>% 
 as.matrix() %>% 
  as.table()

res.ca <- CA(mca, graph = F)

fviz_ca_biplot(res.ca, repel = T)



#MPA

data("poison")

head(poison[,1:7], 3)

pa <- poison[1:55, 5:15]
head(pa[, 1:6], 3)

summary(pa)[, 1:4]
res.mca <- MCA(pa, graph = F)


# RDA for my data

# Load necessary library
library(vegan)

# Load necessary library
library(vegan)

# Example data (replace with your actual data)
# Phytoplankton species data (samples as rows and species as columns)

# RDA analysis with my data
library(vegan)
myrdata <- read_excel("Phyto.xlsx", sheet = "RDA")

tdata <- myrdata %>% 
  select(2:ncol(myrdata)) %>% 
  mutate(across( -pH, ~log10(.+1))) #%>% 
  View()

spdata <- tdata %>% 
  select(1:12) %>% 
  as.matrix()

rownames(spdata) <- c("s1", "s2", "s3", "s4","s5", "s6", "s7", "s8",
                      "w1", "w2", "w3", "w4","w5", "w6", "w7", "w8" )
pdata <-  as.data.frame(spdata)

  
endata <- tdata %>% 
  select(13:ncol(tdata)) %>% 
  as.matrix()


rownames(endata) <- c("s1", "s2", "s3", "s4","s5", "s6", "s7", "s8",
                      "w1", "w2", "w3", "w4","w5", "w6", "w7", "w8" )
edata <- as.data.frame(endata)

seasons <- factor(c("summer", "summer", "winter", "winter"))
rda_result <- rda(pdata ~ ., data = edata)


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



# Network relationship diagram showed the relationship revealing the modular 
#associations among different genus in (a) summer, (b) winter,


library(Hmisc)
library(igraph)
library(ggraph)
library(tidyverse)
library(reshape2)
# Example data - replace this with your actual data
set.seed(123) # For reproducibility
summer_data <- matrix(rnorm(100), ncol=10)
winter_data <- matrix(rnorm(100), ncol=10)
colnames(summer_data) <- paste0("Species", 1:10)
colnames(winter_data) <- paste0("Species", 1:10)

# Calculate Spearman correlations
summer_corr <- rcorr(as.matrix(summer_data), type = "spearman")
winter_corr <- rcorr(as.matrix(winter_data), type = "spearman")

# Extract correlation matrices
summer_corr_matrix <- summer_corr$r
winter_corr_matrix <- winter_corr$r

# Function to convert correlation matrix to edge list and create igraph object
create_network <- function(corr_matrix, threshold = 0.7) {
  corr_matrix[lower.tri(corr_matrix, diag = TRUE)] <- NA  # Remove lower triangle and diagonal
  edges <- melt(corr_matrix, na.rm = TRUE)
  edges <- edges[abs(edges$value) > threshold, ]
  graph <- graph_from_data_frame(edges, directed = FALSE)
  return(graph)
}

# Create networks for summer and winter
summer_network <- create_network(summer_corr_matrix)
winter_network <- create_network(winter_corr_matrix)


# Function to plot network
plot_network <- function(graph, title) {
  ggraph(graph, layout = "fr", weights = abs(E(graph)$value)) + 
    geom_edge_link(aes(edge_alpha = abs(value), edge_color = sign(value)), show.legend = FALSE) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void() +
    ggtitle(title)
}

# Plot for summer
plot_network(summer_network, "Network Relationship Diagram - Summer")

# Plot for winter
plot_network(winter_network, "Network Relationship Diagram - Winter")

# Export networks to Gephi
write_graph(summer_network, file = "summer_network.graphml", format = "graphml")
write_graph(winter_network, file = "winter_network.graphml", format = "graphml")













#working with my actual data

library(readxl)
library(dplyr)
library(tidyverse)
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

#Warer quality parameters
#winter pH
water %>% 
  filter(Parameters == "pH")
set.seed(123)
# Set the parameters for each season and each station
mean_pH_summer <- c(8.1, 8.2, 7.7, 7.9, 8.3, 8.0, 8.4, 8.2)  # Mean pH values for each station in summer
mean_pH_winter <- c(5.9, 7.8, 7.8, 8.1, 7.7, 7.8, 7.3, 7.5)  # Mean pH values for each station in winter

sd_pH_summer <- c(0.11, 0.2, 0.15, 0.1, 0.2, 0.15, 0.05, 0.2)  # Standard deviation of pH values for each station in summer
sd_pH_winter <- c(0.2, 0.15, 0.16, 0.2, 0.15, 0.14, 0.2, 0.15)  # Standard deviation of pH values for each station in winter

num_replicates <- 10  # Number of replicates
num_stations <- length(mean_pH_summer)  # Number of stations

# Generate replicates for each station and season
replicate_data_summer <- matrix(NA, nrow = num_replicates, ncol = num_stations)
replicate_data_winter <- matrix(NA, nrow = num_replicates, ncol = num_stations)

for (i in 1:num_stations) {
  replicate_data_summer[, i] <- rnorm(num_replicates, mean = mean_pH_summer[i], sd = sd_pH_summer[i])
  replicate_data_winter[, i] <- rnorm(num_replicates, mean = mean_pH_winter[i], sd = sd_pH_winter[i])
}

# Convert the matrices to data frames
replicate_df_summer <- as.data.frame(replicate_data_summer)
replicate_df_winter <- as.data.frame(replicate_data_winter)

# Set column names for the stations
colnames(replicate_df_summer) <- paste0("Station_s", 1:num_stations)
colnames(replicate_df_winter) <- paste0("Station_w", 1:num_stations)

pHdata <- cbind(replicate_df_summer, replicate_df_winter)


# DO data
water %>% 
  filter(Parameters == "DO")
set.seed(123)
# Set the parameters for each season and each station
mean_pH_summer <- c(9.13, 8.46, 10.2, 9.78, 10.6, 10.4, 10.21, 10.85)  # Mean pH values for each station in summer
mean_pH_winter <- c(6.89, 6.69, 7.13, 7.11, 6.9, 6.83, 7.01, 7.05)  # Mean pH values for each station in winter

sd_pH_summer <- c(0.31, 0.62, 0.35, 0.41, 0.12, 0.50, 0.25, 0.42)  # Standard deviation of pH values for each station in summer
sd_pH_winter <- c(0.27, 0.35, 0.17, 0.23, 0.35, 0.54, 0.61, 0.43)  # Standard deviation of pH values for each station in winter

num_replicates <- 10  # Number of replicates
num_stations <- length(mean_pH_summer)  # Number of stations

# Generate replicates for each station and season
replicate_data_summer <- matrix(NA, nrow = num_replicates, ncol = num_stations)
replicate_data_winter <- matrix(NA, nrow = num_replicates, ncol = num_stations)

for (i in 1:num_stations) {
  replicate_data_summer[, i] <- rnorm(num_replicates, mean = mean_pH_summer[i], sd = sd_pH_summer[i])
  replicate_data_winter[, i] <- rnorm(num_replicates, mean = mean_pH_winter[i], sd = sd_pH_winter[i])
}

# Convert the matrices to data frames
replicate_df_summer <- as.data.frame(replicate_data_summer)
replicate_df_winter <- as.data.frame(replicate_data_winter)

# Set column names for the stations
colnames(replicate_df_summer) <- paste0("Station_s", 1:num_stations)
colnames(replicate_df_winter) <- paste0("Station_w", 1:num_stations)

DOdata <- cbind(replicate_df_summer, replicate_df_winter)

# Temperature data
water %>% 
  filter(Parameters == "Temp")
set.seed(123)
# Set the parameters for each season and each station
mean_pH_summer <- c(29.4, 28.46, 30.2, 31.78, 30.6, 28.4, 31.29, 28.85)  # Mean pH values for each station in summer
mean_pH_winter <- c(25.89, 26.69, 28.13, 27.11, 26.9, 23.83, 24.01, 26.05)  # Mean pH values for each station in winter

sd_pH_summer <- c(0.51, 0.22, 0.45, 0.11, 0.42, 0.30, 0.45, 0.22)  # Standard deviation of pH values for each station in summer
sd_pH_winter <- c(0.17, 0.25, 0.37, 0.43, 0.15, 0.34, 0.51, 0.33)  # Standard deviation of pH values for each station in winter

num_replicates <- 10  # Number of replicates
num_stations <- length(mean_pH_summer)  # Number of stations

# Generate replicates for each station and season
replicate_data_summer <- matrix(NA, nrow = num_replicates, ncol = num_stations)
replicate_data_winter <- matrix(NA, nrow = num_replicates, ncol = num_stations)

for (i in 1:num_stations) {
  replicate_data_summer[, i] <- rnorm(num_replicates, mean = mean_pH_summer[i], sd = sd_pH_summer[i])
  replicate_data_winter[, i] <- rnorm(num_replicates, mean = mean_pH_winter[i], sd = sd_pH_winter[i])
}

# Convert the matrices to data frames
replicate_df_summer <- as.data.frame(replicate_data_summer)
replicate_df_winter <- as.data.frame(replicate_data_winter)

# Set column names for the stations
colnames(replicate_df_summer) <- paste0("Station_s", 1:num_stations)
colnames(replicate_df_winter) <- paste0("Station_w", 1:num_stations)

Tempdata <- cbind(replicate_df_summer, replicate_df_winter)


# Salinity data
water %>% 
  filter(Parameters == "Salinity")
set.seed(123)
# Set the parameters for each season and each station
mean_pH_summer <- c(26.4, 26.46, 28.2, 28.78, 27.6, 25.4, 27.29, 26.85)  # Mean pH values for each station in summer
mean_pH_winter <- c(31, 33.46, 28.2, 30.78, 29.6, 32.4, 29.29, 32.85)  # Mean pH values for each station in winter

sd_pH_summer <- c(0.61, 0.42, 0.25, 0.31, 0.72, 0.20, 0.55, 0.32)  # Standard deviation of pH values for each station in summer
sd_pH_winter <- c(0.47, 0.35, 0.57, 0.33, 0.55, 0.64, 0.21, 0.43)  # Standard deviation of pH values for each station in winter

num_replicates <- 10  # Number of replicates
num_stations <- length(mean_pH_summer)  # Number of stations

# Generate replicates for each station and season
replicate_data_summer <- matrix(NA, nrow = num_replicates, ncol = num_stations)
replicate_data_winter <- matrix(NA, nrow = num_replicates, ncol = num_stations)

for (i in 1:num_stations) {
  replicate_data_summer[, i] <- rnorm(num_replicates, mean = mean_pH_summer[i], sd = sd_pH_summer[i])
  replicate_data_winter[, i] <- rnorm(num_replicates, mean = mean_pH_winter[i], sd = sd_pH_winter[i])
}

# Convert the matrices to data frames
replicate_df_summer <- as.data.frame(replicate_data_summer)
replicate_df_winter <- as.data.frame(replicate_data_winter)

# Set column names for the stations
colnames(replicate_df_summer) <- paste0("Station_s", 1:num_stations)
colnames(replicate_df_winter) <- paste0("Station_w", 1:num_stations)

Saldata <- cbind(replicate_df_summer, replicate_df_winter)



waterdata <-  cbind(pHdata, DOdata, Tempdata, Saldata)


#write.xlsx(waterdata, "WQdata.xlsx")


wdata <- read_excel("WQdata.xlsx")

# a summary plot with WQ
swq <- wdata %>% 
  pivot_longer(cols = 3:10, names_to = "St", values_to = "values") %>% 
  pivot_wider( names_from = Parameters, values_from = values, values_fn = list(values = mean)) 



write.xlsx(swq, "swq.xlsx")


# for pH

adata <- wdata %>% 
  filter(Parameters == "pH") %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "values")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)



##2.4 plot barplot with errorbars and add cld
ggplot(emmean_cld, aes(Stations, emmean, fill = Seasons, label = .group))+
  geom_bar(stat = "identity", position = "dodge"  )+
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),
                width = 0.2, size = 0.7, position = position_dodge(width = 0.9))+
  geom_text(vjust = -0.7, position = position_dodge(width = 0.9))+
  theme_bw()+
  scale_fill_npg()+
  labs( y = "pH")
# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


# for DO
adata <- wdata %>% 
  filter(Parameters == "DO") %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "values")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)



##2.4 plot barplot with errorbars and add cld
ggplot(emmean_cld, aes(Stations, emmean, fill = Seasons, label = .group))+
  geom_bar(stat = "identity", position = "dodge"  )+
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),
                width = 0.2, size = 0.7, position = position_dodge(width = 0.9))+
  geom_text(vjust = -0.7, position = position_dodge(width = 0.9))+
  theme_bw()+
  scale_fill_npg()+
  labs( y = "DO")
# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)

# for Temp
adata <- wdata %>% 
  filter(Parameters == "Temp") %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "values")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)

##2.4 plot barplot with errorbars and add cld
ggplot(emmean_cld, aes(Stations, emmean, fill = Seasons, label = .group))+
  geom_bar(stat = "identity", position = "dodge"  )+
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),
                width = 0.2, size = 0.7, position = position_dodge(width = 0.9))+
  geom_text(vjust = -0.7, position = position_dodge(width = 0.9))+
  theme_bw()+
  scale_fill_npg()+
  labs( y = "Temperature (C)")
# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


# for salinity
adata <- wdata %>% 
  filter(Parameters == "Salinity") %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "values")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, adata)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)





#2.4 plot barplot with errorbars and add cld
ggplot(emmean_cld, aes(Stations, emmean, fill = Seasons, label = .group))+
  geom_bar(stat = "identity", position = "dodge"  )+
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE),
                width = 0.2, size = 0.7, position = position_dodge(width = 0.9))+
  geom_text(vjust = -0.7, position = position_dodge(width = 0.9))+
  theme_bw()+
  scale_fill_npg()+
  labs( y = "pH")
# on boxplot
adata %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)




#Corelation graph of water quality parameters

library("PerformanceAnalytics")

wdata
coor <- wdata %>% 
  select(3:10)

chart.Correlation(coor, histogram = F, pch= 19)


wdata %>% 
  pivot_longer(cols = 3:10, names_to = "St", values_to = "values") %>% 
  pivot_wider( names_from = Parameters, values_from = values, values_fn = list(values = mean)) %>% 
  select(3:6) %>% 
  chart.Correlation(histogram = T, pch = 10)

wdata %>% 
  filter(Seasons == "Winter") %>% 
  pivot_longer(cols = 3:10, names_to = "St", values_to = "values") %>% 
  pivot_wider( names_from = Parameters, values_from = values, values_fn = list(values = mean)) %>% 
  select(3:6) %>% 
  chart.Correlation(histogram = T, pch = 10)


##############################

# all boxplot
wdata %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "value" ) %>% 
  #filter(Parameters == "Salinity") %>% 
  ggplot(aes(Stations, value, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  facet_wrap(~Parameters, scales = "free_y" )+
  theme_bw()




# all barplot with errorbars

wdata %>% 
  pivot_longer(cols = 3:10, names_to = "Stations", values_to = "value" ) %>% 
  #filter(Parameters == "Salinity") %>% 
  group_by(Seasons, Parameters, Stations) %>% 
  summarise(mean = mean(value),
            sd = sd(value)) %>% 
  ggplot(aes(Stations, mean, fill = Seasons))+
  geom_bar(stat = "identity", position = "dodge", show.legend = F)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
                width = 0.2, size = 0.7, position = position_dodge(width = 0.9))+
  #geom_text(vjust = -0.7, position = position_dodge(width = 0.9))+
  facet_wrap(~Parameters, scales = "free_y" )+
  theme_bw()+
  scale_fill_npg()



# working with phyto abundance and relative abundance data

phyto <- read_excel("Phyto.xlsx", sheet = "phyto")


phyto %>% 
  pivot_longer(cols = 3:14, names_to = "Species", values_to = "Number") %>% 
  group_by(Seasons,Stations) %>% 
  mutate(Number = (Number/ sum(Number)*100)) %>%
  #mutate(Number = sqrt(Number)) %>% 
  ggplot(aes(Stations, Species, size = Number, color = Species))+
  geom_point(show.legend = T)+
  facet_wrap(~Seasons, scales = "free_x")+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 10, family = "Times New Roman"),
        axis.text.y = element_text(face = "italic", ))



phyto %>% 
  pivot_longer(cols = 3:14, names_to = "Species", values_to = "Number" ) %>% 
  ggplot(aes(Stations , Number, fill = Species ))+
  geom_bar(stat = "identity", position = "fill")+
  facet_wrap(~Seasons, scales = "free" )+
  theme_bw()


# working with diversity data

ddata <- read_excel("Phyto.xlsx", sheet = "diversity")

ddata %>% 
  pivot_longer(cols = 4:6, names_to = "Name", values_to = "Values") %>% 
  ggplot(aes(Stations, Values, fill = Seasons))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Summer, scales = "free_y")


# Anova for simpson_1_D diveristy 
simpd <- ddata %>% 
  pivot_longer(cols = 4:6, names_to = "Name", values_to = "values") %>% 
  filter(Summer == "Simpson_1_D")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, simpd)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)

# on boxplot
simpd %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


# Anova for simpson_1_D diveristy 
simpd <- ddata %>% 
  pivot_longer(cols = 4:6, names_to = "Name", values_to = "values") %>% 
  filter(Summer == "Shannon_H")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, simpd)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)

# on boxplot
simpd %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


## Anova for simpson_1_D diveristy 
simpd <- ddata %>% 
  pivot_longer(cols = 4:6, names_to = "Name", values_to = "values") %>% 
  filter(Summer == "Evenness")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, simpd)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)

# on boxplot
simpd %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)


## Anova for Margalef 
simpd <- ddata %>% 
  pivot_longer(cols = 4:6, names_to = "Name", values_to = "values") %>% 
  filter(Summer == "Margalef")
#2.1 anova
mod2 <- aov(values ~ Stations+Seasons, simpd)
summary(mod2)
TukeyHSD(mod2)
#2.2 calculate the summary state
emmean <- emmeans(mod2, specs = c("Stations", "Seasons"))
emmean
#2.3 add cld to emmeans object
emmean_cld <- cld(emmean, Letters = letters)
dim(emmean_cld)

# on boxplot
simpd %>% 
  ggplot(aes(Stations, values, fill = Seasons))+
  geom_boxplot(show.legend = F, outlier.shape = NA)+
  geom_text(data = emmean_cld, aes( Stations, emmean, label = .group),
            position = position_dodge(width = 0.75),
            vjust = -3, hjust = 0.5)








# Generate example time-of-day data (in hours)














